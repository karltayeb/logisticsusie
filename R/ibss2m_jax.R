#' @export
fit_ibss2m_jax <- function(X, y, L = 10,
                           prior_variance = 1.0,
                           estimate_intercept = T,
                           estimate_prior_variance = F,
                           prior_weights = NULL,
                           keep_2m = TRUE) {
  path <- system.file("python", package = "logisticsusie")
  sys <- reticulate::import("sys", convert = FALSE)
  sys$path$append(path)
  reticulate::source_python(system.file("python", "ibss2m.py", package = "logisticsusie"))

  p <- dim(X)[2]
  if (is.null(prior_weights)) {
    prior_weights <- rep(1, p) / p
  }
  tictoc::tic()
  fit <- ibss2m_uvb(X, y,
    L = L,
    prior_variance = prior_variance,
    prior_weights = prior_weights,
    estimate_intercept = estimate_intercept,
    estimate_prior_variance = estimate_prior_variance,
    keep_2m = keep_2m
  )
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}
