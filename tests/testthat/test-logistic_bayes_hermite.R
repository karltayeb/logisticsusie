test_that("compare to adaptive quadrature", {
  sim <- sim_ser()
  x <- sim$X[,1]
  y <- sim$y
  o <- rep(0, length(y))

  expect_equal(logistic_bayes(x, y, o, prior_variance=1)$logZ,
               logistic_bayes_hermite(x, y, o, prior_variance=1, m=33)$logZ)

  expect_equal(logistic_bayes(x, y, o, prior_variance=1)$mu,
               logistic_bayes_hermite(x, y, o, prior_variance=1, m=33)$mu)

  expect_equal(logistic_bayes(x, y, o, prior_variance=1)$var,
               logistic_bayes_hermite(x, y, o, prior_variance=1, m=33)$var)

  expect_equal(logistic_bayes(x, y, o, prior_variance=5)$logZ,
               logistic_bayes_hermite(x, y, o, prior_variance=5, m=33)$logZ)

  expect_equal(logistic_bayes(x, y, o, prior_variance=0.005)$logZ,
               logistic_bayes_hermite(x, y, o, prior_variance=0.005, m=33)$logZ)
})

test_that("compare to laplace", {
  sim <- sim_ser()
  x <- sim$X[,1]
  y <- sim$y
  o <- rep(0, length(y))

  map <- ridge(x, y, o, 1)
  gh <- logistic_bayes_hermite(x, y, o, prior_variance=1, m=1)

  expect_equal(map$mu, gh$mu)
  expect_equal(map$tau, 1/gh$var)
  expect_equal(map$intercept, gh$intercept)
  expect_equal(gh$logZ, map$logZ)
})
