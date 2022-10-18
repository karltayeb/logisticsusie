#' Sigmoid
#' sigmoid function coverts log-odds to probability
#' @param x Log odds
#' @return Returns the probability
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

logSumExp <- matrixStats::logSumExp

softmax <- function(x) {
  return(exp(x - logSumExp(x)))
}

.monotone <- function(v) {
  return(all(tail(v, -1) - head(v, -1) >= 0))
}

#' helper function gets difference in ELBO between iterations
.diff <- function(fit) {
  return(diff(tail(fit$elbo, 2)))
}

#' check if ELBO is converged to some tolerance threshold
.converged <- function(fit, tol = 1e-3) {
  is.converged <- F
  if ((length(fit$elbo) > 1) & .diff(fit) <= tol) {
    is.converged <- T
  }
  return(is.converged)
}


rowCumSum <- function(X) {
  do.call(rbind, apply(X, 1, cumsum, simplify = F))
}


#' Get Credible Sets
#' Wraps susieR::susie_get_cs
#' @export
binsusie_get_cs <- function(fit,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            dedup = TRUE,
                            squared = FALSE,
                            check_symmetric = TRUE,
                            n_purity = 100) {
  res <- list(
    X = fit$data$X,
    alpha = fit$params$alpha,
    V = fit$hypers$prior_variance
  )
  cs <- susieR::susie_get_cs(res,
    X = fit$data$X,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    dedup = dedup,
    squared = squared,
    check_symmetric = check_symmetric,
    n_purity = n_purity
  )
  return(cs)
}


# include the term that cancels out when xi is up-to-date
compute_elbo2 <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi) + 0.5 * omega * (xi^2 - xb2)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# update xi on the fly, so bound is tight
compute_elbo3 <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  xi <- sqrt(xb2)
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}


# implement normal and point mass component distributions

.clamp <- function(v, .upper = 100, .lower = -100) {
  return(pmax(.lower, pmin(v, .upper)))
}





is.odd <- function(x) {
  x %% 2 == 1
}
is.even <- function(x) {
  x %% 2 == 0
}
