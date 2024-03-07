test_that("multiplication works", {
  sim <- sim_ser()
  x <- sim$X[,1]
  y <- sim$y
  fit_fast_glm(x, y, rep(0, length(y)), 0, family='binomial')
})
