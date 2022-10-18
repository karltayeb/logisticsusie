test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("Null model returns 0 PIPs", {
  set.seed(10)
  data <- sim_ser(beta = 0, beta0 = -1)
  fit <- fit.binsusie(data, maxiter = 100, tol = 1e-5)
  fit_pruned <- prune_model(fit, 1)
  fit_pruned <- binsusie_wrapup(fit_pruned, 1e-9)
  expect_equal(fit_pruned$pip, rep(0, length(fit_pruned$pip)))
})
