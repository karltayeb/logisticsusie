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

#' 1d quadrature with a fixed intercept
#' Return the model evidence logp(y| x, b0, b_sigma)
#' which is the marginal likelihood for a logistic model
#' where we integrate over the effect b ~ N(0, b_sigma^2)
fit_quad_1d_fixed_int <- function(x, y, b0, b_sigma, n = 2^10) {
  tictoc::tic()
  res <- quad1d_fixed_intercept2(x, y, b0, 0, b_sigma, n = n)
  tictoc::toc()
  return(res)
}
