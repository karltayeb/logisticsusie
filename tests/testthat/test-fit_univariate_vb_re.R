# these should be identical
test_that("RE model agrees with non-RE model when we dont add a RE", {
  sim <- sim_ser_with_random_effects()
  ser <- with(sim, fit_uvb_ser(X, y, prior_variance = 1.0))
  ser_re0 <- with(sim, fit_uvb_ser_re(X, y, prior_variance = 1.0))
  expect_equal(ser$mu, ser_re0$mu)
  expect_equal(ser$var, ser_re0$var)
  expect_equal(ser$alpha, ser_re0$alpha)
  expect_equal(ser$intercept, ser_re0$intercept)
  expect_equal(ser$lbf, ser_re0$lbf)
  expect_equal(ser$elbo, ser_re0$elbo)
  expect_equal(ser$null_loglik, ser_re0$null_loglik)
  expect_equal(ser$lbf_model, ser_re0$lbf_model)
})


# o <- list(mu=rep(0, 1000), mu2=rep(1, 1000))
# ser1 <- with(sim, fit_uvb_ser_re(X, y, o=o))
#
# o <- add_re(o, ser1, sim$X)
# o$mu
# o$mu2
#
# o <- sub_re(o, ser1, sim$X)
# o$mu
# o$mu2
#
# sim <- sim_susie(length_scale = 10)
# ser1 <- with(sim, fit_uvb_ser_re(X, y))
# plot(ser1$alpha * ser1$mu)
#
# ibss2_fit <- with(sim, ibss2m(X, y, L=3, maxit = 20))
#
# plot(colSums(ibss2_fit$mu * ibss2_fit$alpha))
#
# purrr::map_dbl(ibss2_fit$fits, ~purrr::pluck(.x, 'lbf_model'))
#
# cs_tbl2(ibss2_fit$alpha)
#
# ibss2_fit$mu
