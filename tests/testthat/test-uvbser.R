test_uvbser_N <- function(N = 1) {
  data <- sim_ser(N = N)
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_uvbser(data)
  fit <- update_model(fit, data, fit_prior_variance=T)
  fit <- fit_model(fit, data, fit_prior_variance=T)
  fit2 <- with(data, uvbser(X, y, estimate_prior_variance = F))

  .monotone(fit$elbo)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))

}

# NOTE: this test doesnt actually do anything
# should check that ELBO increases as we e.g. update var0 or include covariates, etc
testthat::test_that("Bernoulli UVB-SER Monotone", {
  for (i in seq(20)) {
    testthat::expect_true(test_uvbser_N(1)$monotone)
  }
})
