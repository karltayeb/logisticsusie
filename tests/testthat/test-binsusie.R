##
# Fit SuSiE tests
##

test_susie_N <- function(N = 1) {
  sim <- sim_susie(N = N, ls=20)
  data <- with(sim, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binsusie(data, L=5)
  fit <- fit_model(fit, data, fit_prior_variance=T, track_elbo=T, max_iter=1000)

  fit2 <- with(sim, binsusie(X, y, N, Z, L=5))

  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SuSiE Monotone", {
  for (i in seq(5)) {
    testthat::expect_true(test_susie_N(1)$monotone)
  }
})


for (i in seq(1, 10)) {
  test_name <- paste("Binomial SuSiE Monotone: N = ", i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_susie_N(i)$monotone)
  })
}

testthat::test_that("p=2 doesn't fail", {
  data <- sim_ser(N = 1)
  data$X <- data$X[, 1:2]
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binsusie(data, 1)
  fit <- fit_model(fit, data, track_elbo=T)
  testthat::expect_true(.monotone(fit$elbo))
})

testthat::test_that("p=1 doesn't fail", {
  data <- sim_ser(N = 1)
  data$X <- data$X[, 1, drop = F]
  data <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binsusie(data, 1)
  fit <- fit_model(fit, data, track_elbo=T)
  testthat::expect_true(.monotone(fit$elbo))
})

testthat::test_that("X a vector throws error", {
  data <- sim_ser(N = 1)
  data$X <- data$X[, 1]
  testthat::expect_error(with(data, binsusie_prep_data(X, y, N, Z)))
})


##
# Fit SuSiE tests
##

test_susie_sparse <- function(N = 1) {
  data <- sim_susie_sparse(N = N)
  data2 <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binsusie(data2, 1)
  fit <- fit_model(fit, data2, track_elbo=T)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SuSiE Monotone (SPARSE)", {
  for (i in seq(5)) {
    testthat::expect_true(test_susie_sparse(1)$monotone)
  }
})

##
# SuSiE with non-zero mu0
##
test_susie_mu0_sigma0 <- function(N = 1, mu0 = 0, var0 = 1, fit_prior_variance = T) {
  set.seed(1)
  data <- sim_susie_sparse(N = N)
  data2 <- with(data, binsusie_prep_data(X, y, N, Z))
  fit <- data_initialize_binsusie(data2, L=5, mu0=mu0, var0=var0)

  fit <- data_initialize_binsusie(data2, L=5, mu0=0, var0=1)
  fit <- fit_model(fit, data2, track_elbo=T)
  .monotone(fit$elbo)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

for (i in seq(1, 10)) {
  test_name <- paste("SuSiE ELBO monotone w/ mu0 != 0")
  testthat::test_that(test_name, {
    testthat::expect_true(test_susie_mu0_sigma0(1, runif(1, -100, 100), 0.5)$monotone)
  })
}
