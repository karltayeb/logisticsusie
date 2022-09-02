testthat::test_that("binomial and logistic bound agree", {
  data <- sim_ser()
  fit <- init.binser(data)
  testthat::expect_equal(jj_bound.binser(fit), jj_bound.logistic(fit))
})
