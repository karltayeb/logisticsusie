#' fit ridge regression
ridge <- function(x, y, o=NULL, prior_variance = 1.){
  # precision = n*lambda, we set the lambda grid so that
  # the largest lambda is our target prior variance
  n <- length(x)

  # Relationship between prior precision and ridge regression penalty in glmnet
  # tau0 = 2 n lambda
  tau0 <- 1/prior_variance
  lambda = 0.5 * rev(exp(seq(log(tau0), log(tau0 * 1e4), length.out = 100))) / n

  # glmnet requres at least two columns
  # we still set intercept = T because otherwise the intercept column gets estimate 0
  # set penalty factor to 0, it sems to give the behavior we want
  X0 <- cbind(rep(1, n), x)
  fit <- glmnet::glmnet(
    X0, y, offset = o,
    intercept=T, alpha=0, standardize = F, penalty.factor = c(0, 1),
    lambda = lambda, family='binomial')


  # summarize MAP
  intercept <- coef(fit, s=lambda[100])[1, 1]
  mu <- coef(fit, s=lambda[100])[3, 1]
  psi <- predict(fit, X0, newoffset=o, s=lambda[100])[,1]
  tau <- sum(sigmoid(psi) * sigmoid(-psi) * x^2) + tau0

  ll <- sum(y * psi - log(1 + exp(psi))) + dnorm(mu, 0, sd=sqrt(1/tau0), log=T)
  ll0 <- -0.5 * fit$nulldev
  #ll <- -0.5 * (fit$nulldev - fit$dev.ratio[100] * fit$nulldev)
  logZ <- ll + 0.5 * log(2 * pi / tau)
  return(list(
    mu = mu, tau=tau, intercept=intercept, prior_variance=prior_variance,
    tau0 = 1/prior_variance,
    ll = ll, ll0=ll0, lr = ll - ll0, logZ = logZ))
}

#' make log likelihood function (as a function of b)
make_log_joint_ridge <- function(x, y, o, ridgefit){
  f <- function(b){
    with(ridgefit, {
      psi <- intercept  + x * b
      if(!is.null(o)){
        psi <- psi + o
      }
      sum((y * psi) - log(1 + exp(psi))) +
        dnorm(b, 0, sd = sqrt(1/tau0), log=T)
    })
  }
  return(Vectorize(f))
}

#' make log likelihood function (as a function of b0, b1)
make_log_joint2d <- function(x, y, o, prior_variance){
  f <- function(b){
    b0 <- b[1]
    b <- b[2]
    psi <- b0  + x * b
    if(!is.null(o)){
      psi <- psi + o
    }
    sum((y * psi) - log(1 + exp(psi))) +
      dnorm(b0, 0, sd = 1000, log=T) +
      dnorm(b, 0, sd = sqrt(prior_variance), log=T)
  }

  f.v <- function(b){
    matrix(apply(b, 2, f), ncol=ncol(b))
  }
  return(f.v)
}

#' make quadratic approx fun (as a function of b)
make_laplace_approx <- function(log_joint, ridgefit){
  g <- function(b){
    with(ridgefit, {
      log_joint(mu) - 0.5 * tau * (b - mu)^2
    })
  }
  return(g)
}

#' create a function f, c -> (f(x) - c)
Shift <- function(f, c){Vectorize(function(x){f(x) - c})}

#' create a function f -> exp(f)
Exp <- function(f){Vectorize(function(x){exp(f(x))})}

#' fit ridge regression and then numerically compute: mean, variance, normalizing constant
#' @export
logistic_bayes <- function(x, y, o=NULL, prior_variance=1., eps=0, width=5){
  # fit ridge regression
  ridgefit <- ridge(x, y, o, prior_variance)

  # f = log p(y, b) -- assuming intercept is fixed at MAP
  f <- make_log_joint_ridge(x, y, o, ridgefit)

  # integrate on b_map +/- width * sqrt(prior_variance)
  # prior_variance > posterior variance
  range <- with(ridgefit, {
    s <- 1/sqrt(tau0)
    c(eps + mu - width*s,  eps + mu + width*s)
  })
  lower <- range[1]
  upper <- range[2]

  if(lower > ridgefit$mu | upper < ridgefit$mu){
    warning('MAP is not in the integration interval')
  }

  #print(c(lower, ridgefit$mu, upper))
  C <- f(ridgefit$mu)
  f2 <- Exp(Shift(f, C))
  logZ <- log(integrate(f2, lower, upper)$value) + C

  D <- exp(C - logZ)
  mu <- integrate(function(b){b * f2(b)}, lower, upper)$value * D
  mu2 <- integrate(function(b){b^2 * f2(b)}, lower, upper)$value * D
  var <- mu2 - mu^2

  return(list(mu = mu, var=var, intercept=ridgefit$intercept, logZ = logZ))
}


#' slow but (hopefully) accurate logistic SER with normal prior
#' @export
logistic_ser <- function(X, y, o=NULL, prior_variance=1.){
  p <- ncol(X)
  n <- length(y)

  if(is.null(o)){
    o <- rep(0, n)
  }
  fit0 <- glm(y ~ 1 + offset(o), family='binomial')
  ll0 <- -0.5 * fit0$null.deviance

  res <- purrr::map(1:p, ~logistic_bayes(
    X[,.x], y, o, prior_variance = prior_variance))

  mu <- purrr::map_dbl(res, ~.x$mu)
  intercept <- purrr::map_dbl(res, ~.x$intercept)
  var <- purrr::map_dbl(res, ~.x$var)
  lbf <- purrr::map_dbl(res, ~.x$logZ) - ll0
  alpha <- exp(lbf - logSumExp(lbf))
  lbf_ser <- logSumExp(lbf) - log(p) # assuming 1/p

  psi <- (X %*% (alpha * mu))[,1]

  return(list(mu=mu, var=var, intercept=intercept, prior_variance=prior_variance, lbf=lbf, lbf_model = lbf_ser, alpha=alpha, psi=psi))
}

