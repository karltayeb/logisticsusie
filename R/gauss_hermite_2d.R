make_quad_points2d <- function(m1, m2){
  quad_points1 <- statmod::gauss.quad(m1, kind='hermite')
  quad_points2 <- statmod::gauss.quad(m2, kind='hermite')
  z1 <- rep(quad_points1$nodes, m2)
  z2 <- rep(quad_points2$nodes, each=m1)
  Z <- t(cbind(z1, z2))

  w <- rep(quad_points1$weights, m2) * rep(quad_points2$weights, each=m1)
  return(list(nodes=Z, weights=w))
}

make_transformed_quad_points2d <- function(H, mu, m1, m2, type='svd', Z=NULL){
  if(is.null(Z)){
    Z <- make_quad_points2d(m1, m2)
  }

  if(type == 'svd'){
    uvt <- svd(solve(-H))
    L <- with(uvt, u %*% diag(sqrt(d)))
    detL <- prod(sqrt(uvt$d))
  } else if(type == 'chol'){
    L <- solve(chol(-H))
    detL <- prod(diag(L))
  } else if(type == 'diag'){
    L <- diag(1/sqrt(diag(-H)))
    detL <- prod(diag(L))
  }

  B <-  (sqrt(2) * L %*% Z$nodes) + mu
  w <- (2 * detL) * Z$weights
  z2 <- colSums(Z$nodes * Z$nodes)
  return(list(nodes=B, weights=w, z2=z2))
}

#' perform 2d gauss hermite quadrature for a given likelihood function
#' B is a list with nodes and weights
.gauss_hermite2d <- function(llfun, B){
  logq <- colSums(B$nodes * B$nodes)
  logZ = logSumExp(llfun(B$nodes[1,], B$nodes[2,]) + B$z2 + log(B$weights))
  return(logZ)
}

#' perform 2d gauss hermite quadrature for a given likelihood function
#' nodes are determined using a local quadratic approximation N(mu, H)
gauss_hermite2d <- function(llfun, H, mu, m1, m2, type='svd', Z=NULL) {
  B <- make_transformed_quad_points2d(H, mu, m1, m2, type=type, Z=Z)
  return(.gauss_hermite2d(llfun, B))
}

compute_hessian <- function(x, y, o, b0, b, prior_variance_b, prior_variance_b0) {
  # compute Hessian
  X0 <- cbind(rep(1, length(x)), x)
  psi <- b0 + x * b + o
  mu <- 1/(1 + exp(-psi))
  f <- sum((y * psi) - log(1 + exp(psi))) +
    dnorm(b, 0, sd = sqrt(prior_variance_b), log=T) +
    dnorm(b0, 0, sd=sqrt(prior_variance_b0), log=T)
  g <- t(X0) %*% (y - mu) - c(b0/prior_variance_b0, b/prior_variance_b)
  H <- -t(X0) %*% (mu * (1 - mu) * X0) - diag(c(1/prior_variance_b0, 1/prior_variance_b))
  return(H)
}

#'
logistic_bayes_hermite2d <- function(x, y, o=NULL, prior_variance=1., prior_variance_b0=1e5, m1=5, m2=5, type='svd', Z=NULL){
  map <- ridge(x, y, o, prior_variance)
  H <- compute_hessian(x, y, o, map$intercept, map$mu, prior_variance, prior_variance_b0)
  llfun <- Vectorize(function(b0, b){
    psi <- b0  + x * b
    if(!is.null(o)){
      psi <- psi + o
    }
    sum((y * psi) - log(1 + exp(psi))) +
      dnorm(b0, 0, sd = sqrt(prior_variance_b0), log=T) +
      dnorm(b, 0, sd = sqrt(prior_variance), log=T)
  })
  logZ <- gauss_hermite2d(llfun, H, c(map$intercept, map$mu), m1, m2, type=type, Z=Z)
  return(list(mu=map$mu, var = -1/H[2,2], intercept=map$intercept, prior_variance=prior_variance, logZ=logZ))
}

#' Bayesian logistic regression with independent gaussian priors on the intercept and effect
#' Numerical integration via Gauss-Hermite quadrature
#' @param X matrix of covariates
#' @param y binary response
#' @param o offset
#' @param prior_variance variance of the prior
#' @param prior_variance_b0 variance of the prior for the intercept
#' @param m1 number of quadrature points for the intercept
#' @param m2 number of quadrature points for the effect
#' @param type type of transformation for the quadrature points, one of 'svd', 'chol', 'diag'
#' @export
logistic_ser_hermite2d <- function(X, y, o=NULL, prior_variance=1., prior_variance_b0=1e5, m1=5, m2=5, type='svd'){
  p <- ncol(X)
  n <- length(y)

  if(is.null(o)){
    o <- rep(0, n)
  }
  fit0 <- glm(y ~ 1 + offset(o), family='binomial')
  ll0 <- -0.5 * fit0$null.deviance
  logZ0 <- ll0 + dnorm(fit0$coefficients[1], 0, sd = sqrt(prior_variance_b0), log=T)

  Z <- make_quad_points2d(m1, m2) # precompute quadrature points
  res <- purrr::map(1:p, ~logistic_bayes_hermite2d(
    X[,.x], y, o, prior_variance = prior_variance, prior_variance_b0 = prior_variance_b0, Z=Z, type=type))

  mu <- purrr::map_dbl(res, ~.x$mu)
  intercept <- purrr::map_dbl(res, ~.x$intercept)
  var <- purrr::map_dbl(res, ~.x$var)
  lbf <- purrr::map_dbl(res, ~.x$logZ) - logZ0
  alpha <- exp(lbf - logSumExp(lbf))
  lbf_ser <- logSumExp(lbf) - log(p) # assuming 1/p

  psi <- (X %*% (alpha * mu))[,1]

  return(list(mu=mu, var=var, intercept=intercept, prior_variance=prior_variance, lbf=lbf, lbf_model = lbf_ser, alpha=alpha, psi=psi))
}
