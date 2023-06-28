test_that("Check dimension of prep_data_for_fastglm", {
  x <- rnorm(100)
  y <- rbinom(100, 1, 0.5)
  o <- rnorm(100)

  # augment and intercept work
  a <- prep_data_for_fastglm(x, y, o, intercept = T, augment = T)
  expect_true(all(dim(a$X) == c(104, 2)))

  # augment works
  b <- prep_data_for_fastglm(x, y, o, intercept = F, augment = T)
  expect_true(all(dim(b$X) == c(104, 1)))

  # no-augment + intercept  works
  c <- prep_data_for_fastglm(x, y, o, intercept = T, augment = F)
  expect_true(all(dim(c$X) == c(100, 2)))

  # no-augment + no-intercept  works
  d <- prep_data_for_fastglm(x, y, o, intercept = F, augment = F)
  expect_true(all(dim(d$X) == c(100, 1)))
})
