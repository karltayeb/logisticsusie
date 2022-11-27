test_that("estimate prior variance improves BFs", {
  sim <- sim_ser()

  # individual lbfs must improve for uvb and glm
  ser <- with(sim, fit_uvb_ser(X, y, estimate_prior_variance = F))
  ser2 <- with(sim, fit_uvb_ser(X, y, estimate_prior_variance = T))
  expect_true(all(ser2$lbf - ser$lbf > 0))

  glm_ser <- with(sim, fit_glm_ser(X, y, prior_variance = 1, estimate_prior_variance = F))
  glm_ser2 <- with(sim, fit_glm_ser(X, y, prior_variance = 1, estimate_prior_variance = T))
  expect_true(all(glm_ser2$lbf - glm_ser$lbf > 0))

  ser_re <- with(sim, fit_uvb_ser_re(X, y, estimate_prior_variance = F))
  ser_re2 <- with(sim, fit_uvb_ser_re(X, y, estimate_prior_variance = T))
  expect_true(all(ser_re2$lbf - ser_re$lbf > 0))

  # model lbf must be improved when prior variance estimate for entire SER
  ser <- with(sim, fit_bin_ser(X, y, estimate_prior_variance = F))
  ser2 <- with(sim, fit_bin_ser(X, y, estimate_prior_variance = T))
  expect_true(all(ser2$lbf_model - ser$lbf_model >= 0))
})
