test_ser_N <- function(N = 1) {
  data <- sim_ser(N = N)
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binser(data)
  fit <- fit_model(fit, data)
  .monotone(fit$elbo)

  compute_psi2(fit, data)
  compute_psi(fit, data)^2

  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone", {
  for (i in seq(20)) {
    testthat::expect_true(test_ser_N(1)$monotone)
  }
})


for (i in seq(1, 10)) {
  test_name <- paste("Binomial SER Monotone: N = ", i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N(i)$monotone)
  })
}


##
# SER With Covariates
##

test_ser_with_covariates_N <- function(N = 1) {
  data <- sim_ser(N = N)
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binser(data)
  fit <- fit_model(fit, data)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER with covariates Monotone", {
  for (i in seq(20)) {
    testthat::expect_true(test_ser_with_covariates_N(1)$monotone)
  }
})


for (i in seq(1, 10)) {
  test_name <- paste("Binomial SER with covariates Monotone: N = ", i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_with_covariates_N(i)$monotone)
  })
}

#' Includes one column of X as a fixed shift
#' So if shift_col = 1, there is basically no effect to estimate
#' If shift_col \neq 1, there are two effects to estimate
test_ser_N_with_shift <- function(N = 1, shift_col = 1, fit_intercept = T, fit_prior_variance = T) {
  data <- sim_ser(N = N)
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  data$shift <- data$X[, shift_col]

  fit <- data_initialize_binser(data)
  fit <- fit_model(fit, data)
  .monotone(fit$elbo)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone w/ shift", {
  set.seed(1)
  for (i in seq(10)) {
    testthat::expect_true(test_ser_N_with_shift(1)$monotone)
  }
})


for (i in seq(1, 10)) {
  test_name <- paste("Binomial SER Monotone w/ shift: N = ", i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N_with_shift(i)$monotone)
  })
}


#' Includes a constant shift
#' basically we want to see that we estimate an intercept that cancels out the shift
test_ser_N_with_constant_shift <- function(N = 1, shift = 1) {
  data <- sim_ser(N = N)
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  data$shift <- shift

  fit <- data_initialize_binser(data)
  fit <- fit_model(fit, data, fit_prior_variance=F)
  .monotone(fit$elbo)

  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone w/ constant shift", {
  for (i in seq(20)) {
    testthat::expect_true(test_ser_N_with_constant_shift(1)$monotone)
  }
})

for (i in seq(1, 10)) {
  test_name <- paste("Binomial SER Monotone w/ constant shift: N = ", i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N_with_constant_shift(i)$monotone)
  })
}


#' Includes a constant random shift (psi2 > psi^2)
#' basically we want to see that we estimate an intercept that cancels out the shift
test_ser_N_with_constant_random_shift <- function(N = 1, shift = 1) {
  data <- sim_ser(N = N)
  data2 <- with(data, binsusie_prep_data(X, y, N, Z))
  data2$shift <- shift
  data2$shift_var <- 5

  fit <- data_initialize_binser(data2)
  fit <- fit_model(fit, data2)

  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone w/ constant random shift", {
  for (i in seq(20)) {
    testthat::expect_true(test_ser_N_with_constant_random_shift(1)$monotone)
  }
})

for (i in seq(1, 10)) {
  test_name <- paste("Binomial SER Monotone w/ constant random shift: N = ", i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N_with_constant_random_shift(i)$monotone)
  })
}

##
# SER with non-zero mu0
##
test_ser_mu0_sigma0 <- function(N = 1, mu0 = 0, var0 = 1, fit_prior_variance = T) {
  set.seed(1)
  data <- sim_ser(N = N)
  data2 <- with(data, binsusie_prep_data(X, y, N, Z))

  fit <- data_initialize_binser(data2, mu0 = mu0, var0 = var0)
  fit <- fit_model(fit, data2, fit_prior_variance=F)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

for (i in seq(1, 10)) {
  test_name <- paste("SER ELBO monotone w/ mu0 != 0")
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_mu0_sigma0(1, runif(1, -100, 100), 0.5)$monotone)
  })
}


