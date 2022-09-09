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
