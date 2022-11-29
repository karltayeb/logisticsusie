test_that("sub_re(add_re(.)) = id", {
  o <- list(mu = rep(0, 100), mu2 = rep(0, 100))
  psi <- rep(0, 100) # rnorm(100)
  ser <- list(psi = psi, psi2 = rgamma(100, 1) + psi^2)
  o2 <- add_re(o, ser, X)
  o3 <- sub_re(o2, ser, X)

  # back to nothing
  all.equal(o$mu, o3$mu)
  all.equal(o$mu2, o3$mu2)

  # add_re positive
  expect_true(all(o2$mu2 - o2$mu^2 > 0))
})

test_that("sub_re(add_re(.)) = id", {
  o <- list(mu = rep(0, 100), mu2 = rep(0, 100))
  psi <- rnorm(100)
  ser <- list(psi = psi, psi2 = rgamma(100, 1) + psi^2)
  o2 <- add_re(o, ser, X)
  o3 <- sub_re(o2, ser, X)

  # back to nothing
  all.equal(o$mu, o3$mu)
  all.equal(o$mu2, o3$mu2)

  # add_re positive
  expect_true(all(o2$mu2 - o2$mu^2 > 0))
})
