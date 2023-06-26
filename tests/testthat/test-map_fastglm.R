test_that('Fast glm looks like glm', {
  sim <- logisticsusie::sim_susie()
  lr1 <- with(sim, map_univariate_regression(X, y))
  lr2 <- with(sim, map_fastglm(X, y))

  expect_equal(lr1$betahat, lr2$betahat)
  expect_equal(lr1$shat2, lr2$shat2)
  expect_equal(lr1$intercept, lr2$intercept)
  expect_equal(lr1$lr, lr2$lr)
})
