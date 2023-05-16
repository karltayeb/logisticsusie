#' @export
update_model.uvbser <- function(fit, data,
                                fit_intercept=TRUE,
                                fit_prior_variance=TRUE,
                                track_elbo = TRUE){

  if(fit_prior_variance){
    fit$var0 <- update_var0.binser(fit)
  }

  newfit <- fit_uvb_ser(data$X, data$y,
              o = data$shift, o2 = (data$shift_var + data$shift^2),
              prior_variance = fit$var0,
              estimate_intercept = fit_intercept,
              prior_weights = fit$pi
  )
  class(newfit) <- c('uvbser', 'binser')
  newfit$var0 <- fit$var0
  newfit$mu0 <- fit$mu0
  newfit$pi <- fit$pi
  newfit$elbo <- c(fit$elbo, newfit$loglik)
  newfit$delta <- sum(fit$alpha * fit$intercept)
  return(newfit)
}

#' @export
fit_model.uvbser <- function(fit, data, ...){
  tictoc::tic()
  fit <- update_model(fit, data, ...)
  timer <- tictoc::toc()
  fit$elapsed_time <- with(timer, toc - tic)
  return(fit)
}

# Initialize -------
#' @export
initialize_uvbser <- function(n, p, p2, mu0=0, var0=1, pi=rep(1/p, p)) {
  # initialize variational parameters
  ser <- list(
    alpha = pi,  # initialize posterior to prior
    mu = rep(mu0, p),
    var = rep(var0, p),
    delta = rep(0, p2),
    pi = pi,
    mu0 = mu0,
    var0 = var0,
    elbo = -Inf
  )
  class(ser) <- c('uvbser', 'binser')
  return(ser)
}

#' @export
data_initialize_uvbser <- function(data, mu0=0, var0=1){
  n <- nrow(data$X)
  p <- ncol(data$X)
  p2 <- ncol(data$Z)
  ser <- initialize_uvbser(n, p, p2, mu0, var0)
  return(ser)
}


#' UVB SER
#'
#' Fit Binomial UVB SER
#' @param X a n x p matrix of covariates
#' @param y an n vector of integer counts, bernoulli/binomial observations
#' @param N the number of binomial trials, defaults to 1, may be a scalar or vector of length n
#' @param Z fixed effect covaraites (including intercept). If null just a n x 1 matrix of ones
#' @param scale if TRUE, scale the columns of X to unit variate
#' @param center if TRUE, center the columns of X to mean zero
#' @param prior_mean the prior mean of each non-zero element of b. Either a scalar or vector of length L.
#' @param prior_variance the prior variance of each non-zero element of b. Either a scalar or vector of length L. If `estimate_prior_variance=TRUE` the value provides an initial estimate of prior variances
#' @param prior_weights prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix
#' @param intercept
#' @param estimate_prior_variance
#' @param returns a fit Binomial SuSiE model, which is compatable with summary functions from `susieR` package
#' @export
uvbser <- function(X,
                   y,
                   o = 0,
                   N = rep(1, length(y)), # number of trials for each
                   Z = NULL,
                   scale = FALSE,
                   center = TRUE,
                   prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
                   prior_variance = 1.0, # TODO: make scaled prior variance?
                   prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect... = hypers$pi
                   intercept = TRUE,
                   estimate_prior_variance = TRUE,
                   check_null_threshold = 0,
                   prior_tol = 1e-09,
                   prune = FALSE,
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   max_iter = 100,
                   tol = 0.001,
                   verbose = FALSE,
                   n_purity = 100) {

  # Prepare data
  data <- binsusie_prep_data(X, y, N, Z, shift=o, scale = scale, center = center)

  # Initialize model
  fit <- data_initialize_uvbser(data, prior_mean, prior_variance)

  # Fit model
  fit <- fit_model(fit, data,
                   fit_intercept=intercept,
                   fit_prior_variance = estimate_prior_variance,
                   track_elbo = TRUE)
  return(fit)
}

#' @export
print.uvbser <- function(fit){
  cs <- get_cs(fit$alpha)
  if(cs$size < 10){
    cs_msg <- paste0('CS = {', paste(cs$cs, collapse=', '), '}')
  } else {
    cs_msg <- paste0('CS = {', paste(head(cs$cs, 10), collapse=', '), ', ... }')
  }

  if(is.null(fit$elbo)){
    fit$elbo <- -Inf
  }
  message <- paste0(
    'UVB SER:',
    '\n ELBO = ', tail(fit$elbo, 1),
    '\n prior_variance = ', fit$var0,
    '\n ', cs_msg
  )
  cat(message)
}
