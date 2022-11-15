test_that("random effect vb agrees with fixed effect vb when there is no offset", {
  sim <- sim_susie(beta = rep(5, 3))

  # models agree when there is no offset
  fit1 <- with(sim, fit_univariate_vb(X[, 1], y))
  fit2 <- with(sim, fit_univariate_vb_re(X[, 1], y))
  fit1$elbos
  expect_equal(fit1$elbos, fit2$elbos)
})


test_that("random effect vb agrees with fixed effect vb when there is a fixed offset", {
  sim <- sim_susie(beta = rep(5, 3))

  # models agree when there is no offset
  fit1 <- with(sim, fit_univariate_vb(X[, 1], y, o = 1))
  fit2 <- with(sim, fit_univariate_vb_re(X[, 1], y, o = list(mu = 1, mu2 = 1)))
  expect_equal(fit1$elbos, fit2$elbos)
})


test_that("uvb_fit_ser = uvb_fit_ser_re NO offset", {
  sim <- sim_susie(beta = rep(5, 3))

  fit1 <- with(sim, fit_uvb_ser(X, y, o = 0))
  fit2 <- with(sim, fit_uvb_ser_re(X, y, o = list(mu = 0, mu2 = 0)))
  expect_equal(fit1$lbf, fit2$lbf)
})


test_that("uvb_fit_ser = uvb_fit_ser_re NO offset", {
  sim <- sim_susie(beta = rep(5, 3))

  fit1 <- with(sim, fit_uvb_ser(X, y, o = 1))
  fit2 <- with(sim, fit_uvb_ser_re(X, y, o = list(mu = 1, mu2 = 1)))
  expect_equal(fit1$lbf, fit2$lbf)
})


# TODO: write test that checks all these agree for all SERs
compute_lbf_model <- function(ser) {
  p <- length(ser$lbf)
  lbf1 <- sum(ser$alpha * ser$lbf) - categorical_kl(ser$alpha, rep(1 / p, p))
  lbf2 <- ser$loglik - ser$null_loglik
  lbf3 <- matrixStats::logSumExp(log(rep(1 / p, p)) + ser$lbf)
  return(c(lbf1, lbf2, lbf3))
}
