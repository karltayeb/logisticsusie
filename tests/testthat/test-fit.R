test_ser_N <- function(N = 1) {
  data <- sim_ser(N = N)
  fit <- fit.binser(data, maxiter = 100, tol = 1e-5)
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
  data <- sim_ser_with_covariates(N = N)
  fit <- fit.binser(data, maxiter = 100, tol = 1e-5)
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
  fit <- init.binser(data)
  shift <- data$X[, shift_col]
  fit <- fit.binser(fit = fit, maxiter = 100, tol = 1e-5, fit_intercept = T, fit_prior_variance = T, shift = shift)
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
  fit <- init.binser(data)
  fit <- fit.binser(fit = fit, maxiter = 100, tol = 1e-5, fit_intercept = T, fit_prior_variance = T, shift = shift)
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


##
# SER with non-zero mu0
##
test_ser_mu0_sigma0 <- function(N = 1, mu0 = 0, sigma0 = 1, fit_prior_variance = T) {
  set.seed(1)
  data <- sim_ser(N = N)
  fit <- init.binser(data)
  fit$hypers$sigma0 <- sigma0
  fit$hypers$mu0 <- mu0

  tol <- 1e-3
  for (i in 1:1000) {
    fit <- iter.binser(fit, fit_intercept = T, fit_prior_variance = T)

    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
    if (.converged(fit, tol)) {
      break
    }
  }
  .monotone(fit$elbo)

  fit$elbo
  fit$params$mu
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

##
# Fit SuSiE tests
##

test_susie_N <- function(N = 1) {
  data <- sim_susie(N = N)
  fit <- fit.binsusie(data, maxiter = 100, tol = 1e-5)
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


##
# Fit SuSiE tests
##

test_susie_sparse <- function(N = 1) {
  data <- sim_susie_sparse(N = N)
  fit <- fit.binsusie(data, maxiter = 100, tol = 1e-5)
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
test_susie_mu0_sigma0 <- function(N = 1, mu0 = 0, sigma0 = 1, fit_prior_variance = T) {
  set.seed(1)
  data <- sim_susie(N = N)
  fit <- init.binsusie(data)
  fit$hypers$sigma0 <- rep(sigma0, length(fit$hypers$sigma0))
  fit$hypers$mu0 <- rep(mu0, length(fit$hypers$mu0))

  tol <- 1e-3
  for (i in 1:1000) {
    fit <- iter.binsusie(fit)

    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.binsusie(fit))
    if (.converged(fit, tol)) {
      break
    }
  }
  .monotone(fit$elbo)

  fit$elbo
  fit$params$mu
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

###
# Two component covariate moderated
###

test_twococomo_N <- function(N = 1) {
  data <- sim_twococomo(N = N)
  fit <- fit.twococomo(data)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}


testthat::test_that("Two CoCoMo Monotone", {
  for (i in seq(5)) {
    testthat::expect_true(test_twococomo_N(1)$monotone)
  }
})


test_twococomo_sparse <- function(N = 1) {
  data <- sim_twococomo_sparse(N = N)
  fit <- fit.twococomo(data)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}


testthat::test_that("Two CoCoMo Monotone (SPARSE)", {
  for (i in seq(5)) {
    testthat::expect_true(test_twococomo_sparse(1)$monotone)
  }
})

###
# Multinomial susie
###

test_mn_susie <- function() {
  data <- sim_mn_susie()
  fit <- fit.mnsusie(data, L = 3, maxiter = 50)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}


testthat::test_that("Multinomial SuSiE", {
  for (i in seq(5)) {
    testthat::expect_true(test_mn_susie()$monotone)
  }
})


###
# More component covariate moderated
###

test_mococomo_N <- function(N = 1) {
  data <- sim_mococomo(N = N)
  fit <- fit.mococomo(data, maxiter = 1000)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}


testthat::test_that("More CoCoMo Monotone", {
  for (i in seq(5)) {
    testthat::expect_true(test_mococomo_N(1)$monotone)
  }
})
