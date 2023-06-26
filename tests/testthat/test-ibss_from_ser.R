test_that("Compare ABF vs corrected ABF in SuSiE IBSS", {
  sim <- logisticsusie::sim_susie()

  ser1 <- with(sim, fit_glm_ser2(X, y, prior_variance = 1, estimate_prior_variance = T))
  ser2 <- with(sim, fit_glm_ser2(X, y, prior_variance = 1, estimate_prior_variance = F))

  ser3 <- with(sim, fit_glm_ser2(X, y, prior_variance = 1, estimate_prior_variance = T, laplace = F))
  ser4 <- with(sim, fit_glm_ser2(X, y, prior_variance = 1, estimate_prior_variance = F, laplace = F))

  ser5 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))

  fit1 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_glm_ser))
  fit2 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_glm_ser2, ))

  fff <- purrr::partial(fit_glm_ser2, estimate_prior_variance = F)
  fit3 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fff))

  fit_null <- with(sim, ibss_from_ser(X, y[sample(1:length(y), replace = F)], 5, prior_variance = 1, ser_function = fff))

  y_null <- sim$y[sample(1:length(sim$y), replace = F)]

  expect_equal(lr1$betahat, lr2$betahat)


  pip1 <- susieR::susie_get_pip(fit1$alpha)
  pip2 <- susieR::susie_get_pip(fit2$alpha)
  fit1$fits[[1]]$lbf - fit2$fits[[1]]$lbf

  fit1 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_uvb_ser))

  fit_uvb_ser2 <- purrr::partial(fit_uvb_ser, mapper=furrr::future_map)
  fit2 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_uvb_ser2))

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


# parallel is slower than sequential here-- probably copying X to each session
fit_parallel <- function(){
  sim <- logisticsusie::sim_susie(p = 5000, n=2000)

  tictoc::tic()
  fit1 <- with(sim, fit_uvb_ser(X, y))
  tictoc::toc()

  fit_uvb_ser2 <- purrr::partial(fit_uvb_ser, mapper=furrr::future_map)
  plan(multisession, workers = 4)

  tictoc::tic()
  fit1 <- with(sim, fit_uvb_ser2(X, y))
  tictoc::toc()
}
