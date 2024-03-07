make_transformed_quad_points <- function(tau, mu, m, Z=NULL){
  if(is.null(Z)){
    Z <- statmod::gauss.quad(m, kind='hermite')
  }
  sigma <- 1/sqrt(tau)
  new_nodes <-  (sqrt(2) * sigma * Z$nodes) + mu
  new_weights <- sqrt(2) * sigma * Z$weights
  z2 <- Z$nodes^2
  return(list(nodes=new_nodes, weights=new_weights, z2=z2))
}

#' perform 2d gauss hermite quadrature for a given likelihood function
#' B is a list with nodes and weights
.gauss_hermite <- function(llfun, B){
  logZ = logSumExp(llfun(B$nodes) + B$z2 + log(B$weights))
  return(logZ)
}

#' perform 2d gauss hermite quadrature for a given likelihood function
#' nodes are determined using a local quadratic approximation N(mu, H)
gauss_hermite <- function(llfun, tau, mu, m, Z=NULL) {
  sigma <- 1/sqrt(tau)
  if(is.null(Z)){
    b <- statmod::gauss.quad.prob(m, dist='normal', mu=mu, sigma=sigma)
  } else {
    b <- list(
      nodes = sigma*Z$nodes + mu,
      weights = Z$weights
    )
  }
  g <- llfun(b$nodes) - dnorm(b$nodes, mu, sigma, log=T)
  logZ = logSumExp(g  + log(b$weights))
  mu <- sum(b$nodes * exp(g - logZ) * b$weights)
  mu2 <- sum(b$nodes^2 * exp(g - logZ) * b$weights)
  return(list(mu=mu, var = mu2 - mu^2, logZ=logZ))
}

#'
logistic_bayes_hermite <- function(x, y, o=NULL, prior_variance=1., m=5, Z=NULL){
  map <- ridge(x, y, o, prior_variance)
  if(m > 1){
    llfun <- Vectorize(function(b){
      psi <- map$intercept  + x * b
      if(!is.null(o)){
        psi <- psi + o
      }
      sum((y * psi) - log(1 + exp(psi))) +
        dnorm(b, 0, sd = sqrt(prior_variance), log=T)
    })
    res <- gauss_hermite(llfun, map$tau, map$mu, m, Z=Z)
    res$intercept <- map$intercept
    res$prior_variance <- prior_variance
  } else{
    res <- map
    res$var <- 1 / res$tau
  }
  return(res)
}

#' Logistic regression with a gaussian prior on the effect
#'   intercept is fixed to the MAP estimate (note: intercept is not regularized, but the effect is)
#'   computed via gauss hermite quadrature
#'   equivalent to Laplace approximation when m=1
#' @param X design matrix
#' @param y binary response
#' @param o offset
#' @param prior_variance prior variance of the effect
#' @param m number of quadrature points
#' @export
logistic_ser_hermite <- function(X, y, o=NULL, prior_variance=1., m=5){
  p <- ncol(X)
  n <- length(y)

  if(is.null(o)){
    o <- rep(0, n)
  }

  fit0 <- glm(y ~ 1 + offset(o), family='binomial')
  ll0 <- -0.5 * fit0$null.deviance
  logZ0 <- ll0

  z <- statmod::gauss.quad.prob(m, dist='normal') # precompute quadrature points
  res <- purrr::map(1:p, ~logistic_bayes_hermite(
    X[,.x], y, o, prior_variance = prior_variance, Z=z))

  mu <- purrr::map_dbl(res, ~.x$mu)
  intercept <- purrr::map_dbl(res, ~.x$intercept)
  var <- purrr::map_dbl(res, ~.x$var)
  lbf <- purrr::map_dbl(res, ~.x$logZ) - logZ0
  alpha <- exp(lbf - logSumExp(lbf))
  lbf_ser <- logSumExp(lbf) - log(p) # assuming 1/p

  psi <- (X %*% (alpha * mu))[,1]

  return(list(mu=mu, var=var, intercept=intercept, prior_variance=prior_variance, lbf=lbf, lbf_model = lbf_ser, alpha=alpha, psi=psi))
}
