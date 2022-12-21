
#' @export
fit_uvb_ser_jax <- function(X, y, o = NULL,
                            prior_variance = 1.0,
                            intercept.init = logodds(mean(y) + 1e-10),
                            estimate_intercept = T,
                            estimate_prior_variance = F,
                            prior_weights = NULL) {
  sys <- reticulate::import("sys")
  sys$path$append(system.file(package = "logisticsusie"))
  reticulate::source_python(system.file("python", "ser.py", package = "logisticsusie"))
  tictoc::tic()
  fit <- fit_uvb_ser(X, y, prior_variance)
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
  reticulate::source_python(system.file("python", "ser.py", package = "logisticsusie"))
  tictoc::tic()
  fit <- fit_tilted_ser(X, y, prior_variance)
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}

#' @export
ibss2m_uvb_jax <- function(X, y, L = 10,
                           prior_variance = 1.0,
                           intercept.init = logodds(mean(y) + 1e-10),
                           estimate_intercept = T,
                           estimate_prior_variance = F,
                           prior_weights = NULL,
                           tol = 1e-5, maxit = 100, keep_2m = TRUE) {
  reticulate::source_python(system.file("python", "ibss2m.py", package = "logisticsusie"))
  tictoc::tic()
  fit <- ibss2m_uvb(X, y, as.integer(L), as.numeric(prior_variance), tol = tol, maxit = as.integer(maxit), keep_2m = as.logical(keep_2m))
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
  reticulate::source_python(system.file("python", "ibss2m.py", package = "logisticsusie"))
  tictoc::tic()
  fit <- ibss2m_tilted(X, y, as.integer(L), as.numeric(prior_variance), tol = tol, maxit = as.integer(maxit), keep_2m = as.logical(keep_2m))
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}
