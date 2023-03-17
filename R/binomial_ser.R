
###
# Single Effect Regression
###

###
# L-indexed getters
###

get_alpha.binser <- function(fit){fit$alpha}
get_pi.binser <- function(fit){fit$pi}
get_mu.binser <- function(fit){fit$mu}
get_var.binser <- function(fit){fit$var}
get_delta.binser <- function(fit){fit$delta}
get_xi.binser <- function(fit){fit$xi}
get_mu0.binser <- function(fit){fit$mu0}
get_var0.binser <- function(fit){fit$var0}

#' Get Binomial SER coefficients
#' Extract regression coefficients from binomial SER fit
#' @param fit Binomial SER object
#' @return Return E[\beta]
coef.binser <- function(fit) {
  b <- Matrix::drop(fit$alpha * fit$mu)
  return(b)
}

#' compute_kl.binser
#' Compute KL divergence between approximate posterior and prior
#' @param fit Binomial SER object
#' @return Return KL(q(\beta) || p(\beta))
compute_kl.binser <- function(fit) {
  return(compute_ser_kl(
    fit$alpha,
    fit$pi,
    fit$mu,
    fit$var,
    fit$mu0,
    fit$var0
  ))
}

#' fixed effects, intercept, shift
compute_Zd.binser <- function(fit, data) {
  Zd <- Matrix::drop(data$Z %*% fit$delta)
  return(Zd)
}

#' compute_Xb.binser
#' Compute expected predicted log odds
#' @param fit Binomial SER object
#' @return Return E[Xb]
compute_Xb.binser <- function(fit, data) {
  Xb <- Matrix::drop(data$X %*% coef(fit))
  return(Xb)
}

compute_psi.binser <- function(fit, data){
  Xb <- compute_Xb(fit, data)
  Zd <- compute_Zd.binser(fit, data)
  psi <- Xb + Zd + data$shift
  return(psi)
}

#' compute_Xb2.binser
#' Compute second moments of predicted log odds
#' @param fit Binomial SER object
#' @return Return E[Xb]
compute_Xb2.binser <- function(fit, data) {
  b2 <- fit$alpha * (fit$mu^2 + fit$var)
  Xb2 <- Matrix::drop(data$X2 %*% b2)
  return(Xb2)
}

compute_psi2.binser <- function(fit, data){
  Xb2 <- compute_Xb2(fit, data)
  Xb <- compute_Xb(fit, data)
  VXb <- Xb2 - Xb^2
  shift <- data$shift
  Zd <- compute_Zd.binser(fit, data)
  psi2 <- (Xb + Zd + data$shift)^2 + VXb + data$shift_var
  return(psi2)
}

#' Compute a scaled posterior mean
#' TODO: add support for non-zero prior mean
#' @param fit a SER object
compute_nu.binser <- function(fit, data) {
  omega <- compute_omega(fit, data)
  kappa <- compute_kappa(data)

  tau0 <- 1. / fit$var0
  prior_mean <- fit$mu0

  Zd <- compute_Zd.binser(fit, data)
  tmp <- kappa - omega * (Zd + data$shift)
  nu <- Matrix::drop(tmp %*% data$X) + (prior_mean * tau0)
  return(nu)
}

#' update intercept and fixed effect covariates
update_delta.binser <- function(fit, data) {
  omega <- compute_omega(fit, data)
  kappa <- compute_kappa(data)
  Xb <- compute_Xb(fit, data) + data$shift
  Z <- data$Z
  delta <- Matrix::drop(Matrix::solve((omega * t(Z)) %*% Z, t(Z) %*% (kappa - omega * Xb)))
  return(delta)
}

#' update SER variational parameters
#' @param nu
#' @param tau
#' @param pi prior probability vector (length p)
update_b_ser <- function(nu, tau, pi) {
  logits <- log(pi) - 0.5 * log(tau) + 0.5 * nu^2 / tau
  logits <- logits - logSumExp(logits)
  alpha <- exp(logits)

  post <- list(
    mu = nu / tau,
    var = 1 / tau,
    alpha = alpha
  )
  return(post)
}

#' update q(b)
update_b.binser <- function(fit, data) {
  nu <- Matrix::drop(compute_nu.binser(fit, data))
  tau0 <- 1 / fit$var0
  tau <- tau0 + fit$tau
  post <- update_b_ser(nu, tau, fit$pi)
  return(post)
}

#' Variational M step for updating prior_variance
#' Optimize ELBO wrt to prior effect standard deviation
update_var0.binser <- function(fit) {
  b2 <- fit$alpha * (fit$mu^2 + fit$var)
  b <- coef(fit)
  mu0 <- fit$mu0
  var0 <- sum(b2 - 2 * b * mu0) + mu0^2
  return(var0)
}


update_model.binser <- function(fit, data,
                                fit_alpha = TRUE,
                                fit_intercept=TRUE,
                                fit_prior_variance=TRUE,
                                fit_xi = TRUE,
                                track_elbo = TRUE){
  # update b
  post <- update_b.binser(fit, data)
  fit$mu <- post$mu
  fit$var <- post$var

  if(fit_alpha){
    fit$alpha <- post$alpha
  }

  # update intercept/fixed effect covariates
  if (fit_intercept) {
    fit$delta <- update_delta.binser(fit, data)
  }
  if (fit_prior_variance) {
    fit$var0 <- update_var0.binser(fit)
  }
  # update xi
  if (fit_xi){
    fit$xi <- update_xi(fit, data)
    fit$tau <- compute_tau(fit, data)
  }

  if(track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }

  return(fit)
}

#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
#' @returns a vector of length N with the lower bound for each data point
compute_jj.binser <- function(fit, data) {
  xi <- fit$xi
  psi <- compute_psi(fit, data)
  kappa <- compute_kappa(data)
  n <- data$N
  bound <- n * log(sigmoid(xi)) + (kappa * psi) - (0.5 * n * xi)
  return(bound)
}

compute_elbo.binser <- function(fit, data) {
  jj <- sum(compute_jj(fit, data)) # E[p(y | w, b)] - KL[q(b) || p(b)]
  kl <- compute_kl(fit) # KL[q(b) || p(b)]
  return(jj - kl)
}


# Initialize -------
initialize_binser <- function(n, p, p2, mu0=0, var0=1, pi=rep(1/p, p)) {
  # initialize variational parameters
  ser <- list(
    alpha = pi,  # initialize posterior to prior
    mu = rep(mu0, p),
    var = rep(var0, p),
    delta = rep(0, p2),
    xi = rep(1e-3, n),
    pi = pi,
    mu0 = mu0,
    var0 = var0,
    tau = rep(1, p), # initialized at end of this function
    elbo = -Inf
  )
  class(ser) <- 'binser'
  return(ser)
}

data_initialize_binser <- function(data, mu0=0, var0=1){
  n <- nrow(data$X)
  p <- ncol(data$X)
  p2 <- ncol(data$Z)
  ser <- initialize_binser(n, p, p2, mu0, var0)
  ser$tau <- compute_tau(ser, data)
  return(ser)
}


#' Binomial SuSiE
#' Fit Binomial SuSiE via coordinate ascent variational inference
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
#' @param s_init a logistic susie object to initialize with, NOTE if non-null, we ignore `prior_mean`, `prior_variance`, and `prior_weights`
#' @param returns a fit Binomial SuSiE model, which is compatable with summary functions from `susieR` package
#' @export
binser <- function(X,
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
                   s_init = NULL, # previous fit with which to initialize NO
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   max_iter = 100,
                   tol = 0.001,
                   verbose = FALSE,
                   n_purity = 100) {

  # Prepare data
  data <- binsusie_prep_data(X, y, N, Z, scale = scale, center = center)

  # Initialize model
  fit <- data_initialize_binser(data, prior_mean, prior_variance)

  # Fit model
  fit <- fit_model(fit, data,
                   fit_intercept=intercept,
                   fit_prior_variance = estimate_prior_variance,
                   track_elbo = TRUE)
  return(fit)
}


print.binser <- function(fit){
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
    'Binomial SER:',
    '\n ELBO = ', tail(fit$elbo, 1),
    '\n prior_variance = ', fit$var0,
    '\n ', cs_msg
  )
  cat(message)
}
