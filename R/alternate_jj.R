explicit_elbo.binser <- function(fit, data) {
  Xb <- compute_Xb(fit, data)
  kappa <- compute_kappa(data)
  xi <- get_xi(fit)
  n <- data$N

  Xb2 <- compute_Xb2(fit)
  omega <- compute_omega(fit, data)

  ll <- -n * log(2) + kappa * Xb - 0.5 * Xb2 * omega
  kl <- pg_kl(n, xi)

  return(ll - kl)
}

#' More explicit, only works for Bernoulli observations
jj_bound.logistic <- function(fit, shift = 0, kidx = NULL) {
  xi <- .get_xi(fit, kidx)
  Xb <- compute_Xb.binser(fit, shift = shift)
  kappa <- compute_kappa(fit, kidx)

  Xb2 <- compute_Xb2.binser(fit, shift = shift)
  omega <- compute_omega(fit)
  bound <- kappa * Xb - 0.5 * xi - 0.25 / xi * tanh(xi / 2) * (Xb2 - xi^2) - log(1 + exp(-xi))
  return(bound)
}

