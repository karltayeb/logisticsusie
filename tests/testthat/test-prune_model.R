test_that("Null model returns 0 PIPs", {
  set.seed(10)
  data <- sim_ser(beta = 0, beta0 = -1)
  fit <- init.binsusie(data)
  fit <- fit.binsusie(fit, maxiter = 100, tol = 1e-5)
  fit_pruned <- prune_model(fit, 1)
  fit_pruned <- binsusie_wrapup(fit_pruned, prior_tol = 1e-3)
  expect_equal(fit_pruned$pip, rep(0, length(fit_pruned$pip)))
})
