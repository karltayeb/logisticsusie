test_that("Compare ABF vs corrected ABF in SuSiE IBSS", {
  sim <- logisticsusie::sim_susie()
  test_gibss <- with(sim, generalized_ibss(X, y, L=5))

  ser1 <- with(sim, fit_glm_ser(X, y, prior_variance = 1, estimate_prior_variance = T))
  ser2 <- with(sim, fit_glm_ser(X, y, prior_variance = 1, estimate_prior_variance = F))

  ser3 <- with(sim, fit_glm_ser(X, y, prior_variance = 1, estimate_prior_variance = T, laplace = F))
  ser4 <- with(sim, fit_glm_ser(X, y, prior_variance = 1, estimate_prior_variance = F, laplace = F))

  ser5 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))

  fit1 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_glm_ser))
  fit2 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = fit_glm_ser, ))

  fff <- purrr::partial(fit_glm_ser, estimate_prior_variance = F)
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

test_that("Generalized IBSS (binomial)  works with probabilities", {
  sim <- logisticsusie::sim_susie(beta0 = -1)
  y2 <- pmax(pmin(sim$y, 0.99), 0.01)
  test_gibss <- with(sim, generalized_ibss(X, y2, L=5))
  psi <- predict(test_gibss, sim$X)
})

test_that("Init works", {
  sim <- logisticsusie::sim_ser()
  fit <- with(sim, ibss_from_ser(X, y, prior_variance = 1, ser_function = fit_glm_ser, maxit=20))

  # suppress "maximum iterations reached" warning
  suppressWarnings({
    fit1 <- with(sim, ibss_from_ser(X, y, prior_variance = 1, ser_function = fit_glm_ser, maxit=1))
    fit_from_init <- with(sim, ibss_from_ser(X, y, prior_variance = 1, ser_function = fit_glm_ser, init=fit1, maxit=1))
    fit2 <- with(sim, ibss_from_ser(X, y, prior_variance = 1, ser_function = fit_glm_ser, maxit=3))
  })

  expect_equal(fit2$alpha, fit_from_init$alpha)
  expect_equal(fit2$mu, fit_from_init$mu)
  expect_equal(fit2$var, fit_from_init$var)
})


test_that("Compare ABF vs correct ABF in SER", {
  sim <- logisticsusie::sim_ser()
  fit1 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))
  fit2 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))
})


test_that("Compare ABF vs correct ABF in SER", {
  sim <- logisticsusie::sim_ser()
  fit1 <- with(sim, binser(X, y))
  fit1 <- with(sim, fit_glm_ser(X, y, prior_variance = 1))
  fit2 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = binser))
  fit3 <- with(sim, ibss_from_ser(X, y, 5, prior_variance = 1, ser_function = uvbser))
})

test_that("ibss_from_ser gives same result with num_cores > 1",{
  library(survival)
  data(survival_demo)
  X <- survival_demo$geno
  pheno <- with(survival_demo,Surv(time,status))
  set.seed(1)
  fit1 <- with(survival_demo,
               ibss_from_ser(X,pheno,L = 3,
                             ser_function = ser_from_univariate(surv_uni_fun),
                             num_cores = 1))
  set.seed(1)
  fit2 <- with(survival_demo,
               ibss_from_ser(X,pheno,L = 3,
                             ser_function = ser_from_univariate(surv_uni_fun),
                             num_cores = 4))
  fit1$elapsed_time <- 0
  fit2$elapsed_time <- 0
  expect_equal(fit1,fit2)
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
