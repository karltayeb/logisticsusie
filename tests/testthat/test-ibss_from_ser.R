test_that("Compare ABF vs corrected ABF in SuSiE IBSS", {
  sim <- logisticsusie::sim_susie()
  fit1 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_glm_ser))
  fit2 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_glm_ser2))

  pip1 <- susieR::susie_get_pip(fit1$alpha)
  pip2 <- susieR::susie_get_pip(fit2$alpha)
  fit1$fits[[1]]$lbf - fit2$fits[[1]]$lbf
})


test_that("Compare ABF vs correct ABF in SER", {
  sim <- logisticsusie::sim_ser()
  fit1 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))
  fit2 <- with(sim, fit_glm_ser2(X, y, prior_variance = 1))
})


test_that("Compare ABF vs correct ABF in SER", {
  sim <- logisticsusie::sim_ser()
  fit1 <- with(sim, binser(X, y))
  fit1 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))
  fit2 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = binser))
  fit3 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = uvbser))
})
