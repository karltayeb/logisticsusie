

fit_uvb_ser_re_jax <- function(X, y, o = NULL,
                               prior_variance = 1.0,
                               intercept.init = logodds(mean(y) + 1e-10),
                               estimate_intercept = T,
                               estimate_prior_variance = F,
                               prior_weights = NULL) {
  reticulate::source_python("~/R/logisticsusie/inst/python/univariate_vb.py")
  fit <- fit_uvb_ser_jax2(X, y, prior_variance)
  return(fit)
}

ibss2m_jax <- function(X, y, L = 10,
                       prior_variance = 1.0,
                       intercept.init = logodds(mean(y) + 1e-10),
                       estimate_intercept = T,
                       estimate_prior_variance = F,
                       prior_weights = NULL, tol = 1e-5) {
  reticulate::source_python("~/R/logisticsusie/inst/python/ibss2m.py")
  fit <- ibss2m_jax(X, y, as.integer(L), as.numeric(prior_variance), tol = tol)
  return(fit)
}
