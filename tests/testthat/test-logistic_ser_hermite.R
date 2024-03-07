test_that("Check similarity of hermite a adaptive schemes", {
  sim <- sim_ser()
  tictoc::tic()
  a <- with(sim, logistic_ser_hermite(X, y, m=11, prior_variance = 1))
  tictoc::toc()

  tictoc::tic()
  b <- with(sim, logistic_ser(X, y, prior_variance = 1))
  tictoc::toc()

  par(mfrow=c(1,2))
  plot(a$lbf, b$lbf); abline(0, 1)
  plot(a$mu, b$mu);abline(0, 1)
})

test_that("Check similarity of hermite a adaptive schemes", {
  sim <- sim_ser()
  tictoc::tic()
  a <- with(sim, logistic_ser_hermite2d(X, y, m1=1, m2=1, type='diag', prior_variance = 1))
  tictoc::toc()

  tictoc::tic()
  b <- with(sim, logistic_ser(X, y, prior_variance = 1))
  tictoc::toc()

  par(mfrow=c(1,2))
  plot(a$lbf, b$lbf); abline(0, 1)
  plot(a$mu, b$mu);abline(0, 1)
})


