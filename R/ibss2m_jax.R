#' @export
ibss2m_uvb_jax <- function(X, y, L = 10,
                           prior_variance = 1.0,
                           intercept.init = logodds(mean(y) + 1e-10),
                           estimate_intercept = T,
                           estimate_prior_variance = F,
                           prior_weights = NULL,
                           tol = 1e-5, maxit = 100, keep_2m = TRUE) {
  path <- system.file("python", package = "logisticsusie")
  sys <- reticulate::import("sys", convert = FALSE)
  sys$path$append(path)
  reticulate::source_python(system.file("python", "ibss2m.py", package = "logisticsusie"))

  tictoc::tic()
  fit <- ibss2m_uvb(X, y, as.integer(L), as.numeric(prior_variance), tol = tol, maxit = as.integer(maxit), keep_2m = as.logical(keep_2m))
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}
