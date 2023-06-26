test_that("fit_glm_ser2 matches reference implementation", {
  sim <- logisticsusie::sim_susie()

  ser1 <- with(sim, fit_glm_ser(X, y))
  ser2 <- with(sim, fit_glm_ser2(X, y, laplace = F, estimate_prior_variance = F))
  ser3 <- with(sim, fit_glm_ser2(X, y, laplace = T, estimate_prior_variance = F))

  expect_equal(ser1$alpha, ser2$alpha)
  expect_equal(ser1$mu, ser2$mu)
  expect_equal(ser1$var, ser2$var)
  expect_equal(ser1$lbf, ser2$lbf)
  expect_equal(ser1$lbf_model, ser2$lbf_model)
  expect_equal(ser1$prior_variance, ser2$prior_variance)
  expect_equal(ser1$betahat, ser2$betahat)
  expect_equal(ser1$shat2, ser2$shat2)
  expect_equal(ser1$intercept, ser2$intercept)
})


test_that("fit_glm_ser2 returns sensible log Laplace ABFs", {
  set.seed(1)
  sim <- logisticsusie::sim_ser()

  laplace <- purrr::partial(fit_glm_ser2, laplace=T, estimate_prior_variance=F)

  # when we fit the SER there is strong evidence for an effect
  ser_laplace <- with(sim, laplace(X, y))
  expect_true(ser_laplace$lbf_model > 100)

  # when we include the explanatory variable as a fixed offset
  # there is no advantave to adding an effect
  ser_laplace_offset <- with(sim, laplace(X, y, o = X[, 1]))
  expect_true(ser_laplace_offset$lbf_model < 0)

  # there is only one effect, which will get picked up by the first ser
  # the rest should have a negative log BF for this simulation setting
  # fit_laplace <- with(sim, ibss_from_ser(X, y, L=5, ser_function = laplace, maxit = 1))
  # expect_true(all(tail(purrr::map_dbl(fit_laplace$fits, ~.$lbf_model),4) < 0))
})


test_that("fastglm is fast?", {
  set.seed(1)
  sim <- logisticsusie::sim_ser()

  slow <- purrr::partial(fit_glm_ser2, laplace=T, estimate_prior_variance=F, glm_mapper=map_univariate_regression)
  fast <- purrr::partial(fit_glm_ser2, laplace=T, estimate_prior_variance=F, glm_mapper=map_fastglm)

  slow_fit <- with(sim, slow(X, y))
  fast_fit <- with(sim, fast(X, y))

  slow_bench <- microbenchmark::microbenchmark(with(sim, slow(X, y)))
  fast_bench <- microbenchmark::microbenchmark(with(sim, fast(X, y)))

  # when we fit the SER there is strong evidence for an effect
  ser_laplace <- with(sim, laplace(X, y))
  expect_true(ser_laplace$lbf_model > 100)

  # when we include the explanatory variable as a fixed offset
  # there is no advantave to adding an effect
  ser_laplace_offset <- with(sim, laplace(X, y, o = X[, 1]))
  expect_true(ser_laplace_offset$lbf_model < 0)

  # there is only one effect, which will get picked up by the first ser
  # the rest should have a negative log BF for this simulation setting
  # fit_laplace <- with(sim, ibss_from_ser(X, y, L=5, ser_function = laplace, maxit = 1))
  # expect_true(all(tail(purrr::map_dbl(fit_laplace$fits, ~.$lbf_model),4) < 0))
})

