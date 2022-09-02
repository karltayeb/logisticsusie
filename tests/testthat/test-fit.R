test_ser_N <- function(N=1){
  data <- sim_ser(N = N)
  fit <- fit.binser(data, maxiter = 100, tol=1e-5)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone", {
  for(i in seq(20)){
    testthat::expect_true(test_ser_N(1)$monotone)
  }
})


for(i in seq(1, 10)){
  test_name <- paste('Binomial SER Monotone: N = ', i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N(i)$monotone)
  })
}


##
# SER With Covariates
##

test_ser_with_covariates_N <- function(N=1){
  data <- sim_ser_with_covariates(N = N)
  fit <- fit.binser(data, maxiter = 100, tol=1e-5)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER with covariates Monotone", {
  for(i in seq(20)){
    testthat::expect_true(test_ser_with_covariates_N(1)$monotone)
  }
})


for(i in seq(1, 10)){
  test_name <- paste('Binomial SER with covariates Monotone: N = ', i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_with_covariates_N(i)$monotone)
  })
}

#' Includes one column of X as a fixed shift
#' So if shift_col = 1, there is basically no effect to estimate
#' If shift_col \neq 1, there are two effects to estimate
test_ser_N_with_shift <- function(N=1, shift_col=1, fit_intercept=T, fit_prior_variance=T){
  data <- sim_ser(N = N)
  fit <- init.binser(data)
  fit$hypers$shift <- data$X[,shift_col]

  fit <- fit.binser(fit=fit, maxiter = 100, tol=1e-5, fit_intercept = T, fit_prior_variance = T)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone w/ shift", {
  set.seed(1)
  for(i in seq(10)){
    testthat::expect_true(test_ser_N_with_shift(1)$monotone)
  }
})


for(i in seq(1, 10)){
  test_name <- paste('Binomial SER Monotone w/ shift: N = ', i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N_with_shift(i)$monotone)
  })
}



#' Includes a constant shift
#' basically we want to see that we estimate an intercept that cancels out the shift
test_ser_N_with_constant_shift <- function(N=1, shift=1){
  data <- sim_ser(N = N)
  fit <- init.binser(data)
  fit$hypers$shift <- shift
  fit <- fit.binser(fit=fit, maxiter = 100, tol=1e-5, fit_intercept = T, fit_prior_variance = T)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

testthat::test_that("Bernoulli SER Monotone w/ constant shift", {
  for(i in seq(20)){
    testthat::expect_true(test_ser_N_with_constant_shift(1)$monotone)
  }
})

for(i in seq(1, 10)){
  test_name <- paste('Binomial SER Monotone w/ constant shift: N = ', i)
  testthat::test_that(test_name, {
    testthat::expect_true(test_ser_N_with_constant_shift(i)$monotone)
  })
}




# DEbug
if(F){
  # failing case
  set.seed(3)
  test_ser_N_with_shift()$fit$elbo

  # isolate
  set.seed(3)
  data <- sim_ser(N = N)
  fit <- init.binser(data)
  fit$hypers$shift <- data$X[,shift_col]
  fit <- fit.binser(fit=fit, maxiter = 2, tol=1e-5, fit_intercept = F, fit_prior_variance = F)
  fit$elbo

  post <- update_b.binser(fit)

  old_mu <- fit$params$mu
  old_var <- fit$params$var
  old_alpha <- fit$params$alpha
  old_jj <- jj_bound.binser(fit)


  new_mu <- post$mu
  new_var <- post$var
  new_alpha <- post$alpha

  fit$params$mu <- post$mu
  fit$params$var <- post$var

  new_jj <- jj_bound.binser(fit)

  sum(old_alpha * normal_kl(
    old_mu, old_var, 0, 1))

  sum(old_alpha * normal_kl(
    new_mu, new_var, 0, 1))


  fit$params$mu <- post$mu
  fit$params$var <- post$var
  fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
  fit$elbo

  fit$params$alpha <- post$alpha
  fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
  fit$elbo

  # update xi
  fit$params$xi <- update_xi.binser(fit)
  fit$params$tau <- compute_tau(fit)
  fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
  fit$elbo
}
