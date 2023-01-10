#' @export
fit_uvb_ser_jax <- function(X, y, o = NULL,
                            prior_variance = 1.0,
                            intercept.init = logodds(mean(y) + 1e-10),
                            estimate_intercept = T,
                            estimate_prior_variance = F,
                            prior_weights = NULL) {
  if (is.null(o)) {
    o <- 0
  }

  path <- system.file("python", package = "logisticsusie")
  sys <- reticulate::import("sys", convert = FALSE)
  sys$path$append(path)
  reticulate::source_python(system.file("python", "ser.py", package = "logisticsusie"))

  tictoc::tic()
  fit <- fit_uvb_ser2(X, y, prior_variance, o = o, o2 = 0)
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}

#' @export
fit_tilted_ser_jax <- function(X, y, o = NULL,
                               prior_variance = 1.0,
                               intercept.init = logodds(mean(y) + 1e-10),
                               estimate_intercept = T,
                               estimate_prior_variance = F,
                               prior_weights = NULL) {
  path <- system.file("python", package = "logisticsusie")
  sys <- reticulate::import("sys", convert = FALSE)
  sys$path$append(path)
  reticulate::source_python(system.file("python", "ser.py", package = "logisticsusie"))

  tictoc::tic()
  fit <- fit_tilted_ser2(X, y, prior_variance)
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}


#' @export
ibss2m_tilted <- function(X, y, L = 10,
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
  fit <- ibss2m_tilted(X, y, as.integer(L), as.numeric(prior_variance), tol = tol, maxit = as.integer(maxit), keep_2m = as.logical(keep_2m))
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}
