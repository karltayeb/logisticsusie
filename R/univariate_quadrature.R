# Compute likelihood conditional on effect size y ~ Bernoulli(sigmoid(x*b + delta))
conditional_likelihood <- function(x, y, b, delta = 0, log = F) {
  p <- sigmoid(x * b + delta)
  ll <- sum(dbinom(y, 1, p, log = T))
  if (!log) {
    ll <- exp(ll) # rescale for numerical stability?
  }
  return(ll)
}

# compute quadrature integration over gaussian
quad1d_fixed_intercept <- function(x, y, b0, b_mu, b_sigma, n = 128) {
  # inner integral over b
  quad_b <- function(b) {
    res <- conditional_likelihood(x, y, b, delta = b0, log = T)
    return(res)
  }

  # integrate over b0
  q <- statmod::gauss.quad.prob(n = n, dist = "normal", mu = b_mu, sigma = b_sigma)
  integrand <- purrr::map_dbl(q$nodes, quad_b)
  res <- matrixStats::logSumExp(integrand + log(q$weights))
  return(res)
}

# compute quadrature integration over uniform
quad1d_fixed_intercept2 <- function(x, y, b0, b_mu, b_sigma, n = 128, q = NULL) {
  # inner integral over b
  quad_b <- function(b) {
    res <- conditional_likelihood(x, y, b, b0, log = T) + dnorm(b, b_mu, b_sigma, log = T)
    return(res)
  }

  # integrate over b0
  q <- statmod::gauss.quad.prob(n = n, l = b_mu - 10 * b_sigma, u = b_mu + 10 * b_sigma)
  integrand <- purrr::map_dbl(q$nodes, quad_b)
  res <- matrixStats::logSumExp(integrand + log(q$weights)) + log(20 * b_sigma)
  return(res)
}

#' Gaussian quadrature
#'
#' approximate log(\int f(x)dx) by quadrature on \int f(x)/p(x) p(x)
#' where p(x) is a normal density, we approximate this integral
#' with a gaussian quadrature rule via `statmod::gauss.quad.prob`
#' done on a log scale for numerical stability (f is positive)
#' the goal is to choose p(x) such that f(x)/p(x) is not too variable near the mode
#' @param mu mean of gaussian quadrature
#' @param sd sd of gaussian quadrature
gauss_quad2 <- function(log_f, mu, sd, n=32){
  q <- statmod::gauss.quad.prob(n, dist='normal', mu=mu, sigma=sd)
  log_integrand <- function(x){log_f(x) - dnorm(x, mu, sd, log=T)}
  return(matrixStats::logSumExp(log_integrand(q$nodes) + log(q$weights)))
}

# log p(y | x, beta0, beta) as a function of beta
make_log_joint <- function(x, y, beta0, sigma2){
  ll_base <- function(beta){
    psi <- beta0 + x * beta
    ll <- sum(psi*y - log(1 + exp(psi))) + dnorm(beta, sd = sqrt(sigma2), log=T)
    return(ll)
  }
  log_joint <- Vectorize(ll_base)
  return(log_joint)
}

#' 1d quadrature with a fixed intercept
#' Return the model evidence logp(y| x, b0, b_sigma)
#' which is the marginal likelihood for a logistic model
#' where we integrate over the effect b ~ N(0, b_sigma^2)
compute_log_marginal <- function(x, y, b0, mu, var, prior_variance, n = 32, verbose=F) {
  if(verbose){
    message(glue::glue('compute log marginal via quadrature n = {n}'))
  }
  log_joint <- make_log_joint(x, y, b0, prior_variance)
  res <- gauss_quad2(log_joint, mu, sqrt(var), n=n)
  return(res)
}


#' Adaptive quadrature
#' Double the number of quadrature points until change is small
compute_log_marginal_adaptive <- function(x, y, b0, mu, var, prior_variance, n_init = 8, tol=1e-3, verbose=T) {
  log_joint <- make_log_joint(x, y, b0, prior_variance)
  n <- n_init
  log_marginal <- -Inf
  while(T){
    if(verbose){
      message(glue::glue('Computing quadrature rule, n = {n}'))
    }
    log_marginal_new <- gauss_quad2(log_joint, mu, sqrt(var), n=n)
    diff <- abs(log_marginal_new - log_marginal)
    log_marginal <- log_marginal_new
    n <- n * 2
    if(diff < tol){
      break
    }
  }
  return(log_marginal)
}

#' Fit quad SER
#'
#' posterior mean and variance are approximated using an asymptotic normal approximation
#'  the same as in `fit_glm_ser2`
#' the BFs are computed via quadrature, with a quadrature rule
#'  centered on the approximate MAP
#' @param X n x p matrix of covariates
#' @param y n vector of responses
#' @param o n vector of fixed offsets
#' @param glm_ser_args arguments to be passed to `fit_glm_ser2`
#' @param glm_ser provide a fit glm_ser object to avoid refitting
#' @param n number of quadrature points
#' @param verbose print messages during fitting
#' @export
fit_quad_ser <- function(X, y, o = NULL, glm_ser_args = list(), glm_ser=NULL, n=16, verbose=F) {
  p <- dim(X)[2]
  if(is.null(glm_ser)){
    glm_ser_args = c(list(X=X, y=y, o=o), glm_ser_args)
    glm_ser <- rlang::exec(fit_glm_ser2, !!!glm_ser_args)
  }

  prior_variance = glm_ser$prior_variance
  null_loglik <- sum(dbinom(y, 1, mean(y), log = T))
  lbf <- purrr::map_dbl(1:p, ~ compute_log_marginal(
    x = X[, .x],
    y = y,
    b0 = glm_ser$intercept[.x],
    mu = glm_ser$mu[.x],
    var = glm_ser$var[.x],
    prior_variance = prior_variance,
    n = n,
    verbose = verbose
  )) - null_loglik

  alpha <- exp(lbf - logSumExp(lbf))
  lbf_model <- sum(alpha * lbf) - categorical_kl(alpha, rep(1 / p, p))
  loglik <- lbf_model + null_loglik

  quad_ser <- glm_ser
  quad_ser$alpha <- alpha
  quad_ser$lbf <- lbf
  quad_ser$lbf_model <- lbf_model
  class(quad_ser) <- c('quad_ser', 'asymptotic_ser')
  return(quad_ser)
}

