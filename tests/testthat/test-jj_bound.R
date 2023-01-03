testthat::test_that("binomial and logistic bound agree", {
  data <- sim_ser()
  fit <- fit.binser(data)
  testthat::expect_equal(jj_bound.binser(fit), jj_bound.logistic(fit))
})

testthat::test_that("Explicit ELBO and  JJ bound agree", {
  data <- sim_ser()
  fit <- fit.binser(data)

  a <- jj_bound.binser(fit)
  b <- explicit_elbo.binser(fit)
  testthat::expect_equal(a, b)
})
