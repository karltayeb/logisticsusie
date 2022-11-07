
###
# Single Effect Regression
###

###
# L-indexed getters
###

#' idx is for indexing when we run susie
.get_alpha <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    alpha <- fit$params$alpha
  } else {
    alpha <- fit$params$alpha[idx, ]
  }
  return(alpha)
}

.get_mu <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    mu <- fit$params$mu
  } else {
    mu <- fit$params$mu[idx, ]
  }
  return(mu)
}

.get_var <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    var <- fit$params$var
  } else {
    var <- fit$params$var[idx, ]
  }
  return(var)
}

.get_pi <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    pi <- fit$hypers$pi
  } else {
    pi <- fit$hypers$pi[idx, ]
  }
  return(pi)
}

.get_delta <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    delta <- fit$params$delta
  } else {
    delta <- fit$params$delta[idx, ]
  }
  return(delta)
}

.get_prior_mean <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    prior_mean <- fit$hypers$prior_mean
  } else {
    prior_mean <- fit$hypers$prior_mean[idx]
  }
  return(prior_mean)
}

.get_var0 <- function(fit, idx = NULL) {
  if (is.null(idx)) {
    var0 <- fit$hypers$prior_variance
  } else {
    var0 <- fit$hypers$prior_variance[idx]
  }
  return(var0)
}

###
# K-indexed getters
###

.get_y <- function(fit, kidx = NULL) {
  if (is.null(fit$kidx)) {
    y <- fit$data$y
  } else {
    y <- fit$data$y[, fit$kidx]
  }
  return(y)
}

.get_N <- function(fit, kidx = NULL) {
  if (is.null(fit$kidx)) {
    N <- fit$data$N
  } else {
    N <- fit$data$N[, fit$kidx]
  }
  return(N)
}

.get_xi <- function(fit, kidx = NULL) {
  # if(is.null(fit$kidx)){
  #   xi <- fit$params$xi
  # } else{
  #   xi <- fit$params$xi[, fit$kidx]
  # }
  xi <- fit$params$xi
  return(xi)
}


#' Get Binomial SER coefficients
#' Extract regression coefficients from binomial SER fit
#' @param fit Binomial SER object
#' @return Return E[\beta]
coef.binser <- function(fit, idx = NULL) {
  b <- Matrix::drop(.get_alpha(fit, idx) * .get_mu(fit, idx))
  return(b)
}

#' KL(q(b, \omega) || p(\beta, \omega))
.compute_ser_kl <- function(alpha, pi, mu, var, prior_mean, var0) {
  kl <- categorical_kl(alpha, pi)
  kl <- kl + sum(alpha * normal_kl(
    mu, var, prior_mean, var0
  ))
  return(kl)
}

#' compute_kl.binser
#' Compute KL divergence between approximate posterior and prior
#' @param fit Binomial SER object
#' @return Return KL(q(\beta) || p(\beta))
compute_kl.binser <- function(fit) {
  return(.compute_ser_kl(
    .get_alpha(fit),
    .get_pi(fit),
    .get_mu(fit),
    .get_var(fit),
    .get_prior_mean(fit),
    .get_var0(fit)
  ))
}

#' fixed effects, intercept, shift
compute_Zd.binser <- function(fit, idx = NULL, shift = 0) {
  Zd <- Matrix::drop(fit$data$Z %*% .get_delta(fit, idx)) + shift
  return(Zd)
}

#' compute_Xb.binser
#' Compute expected predicted log odds
#' @param fit Binomial SER object
#' @return Return E[Xb]
compute_Xb.binser <- function(fit, idx = NULL, shift = 0) {
  Xb <- Matrix::drop(fit$data$X %*% coef.binser(fit, idx))
  Zd <- compute_Zd.binser(fit, idx, shift)
  Xb <- Xb + Zd
  return(Xb)
}

#' compute_Xb2.binser
#' Compute second moments of predicted log odds
#' @param fit Binomial SER object
#' @return Return E[Xb]
compute_Xb2.binser <- function(fit, idx = NULL, shift = 0) {
  Xb <- Matrix::drop(fit$data$X %*% coef.binser(fit, idx))
  Zd <- compute_Zd.binser(fit, idx, shift)

  b2 <- .get_alpha(fit, idx) * (.get_mu(fit, idx)^2 + .get_var(fit, idx))
  Xb2 <- Matrix::drop(fit$data$X2 %*% b2)

  Xb2 <- Xb2 + 2 * Xb * Zd + Zd^2
  return(Xb2)
}

#' Compute E[y - N/2]
compute_kappa <- function(fit, kidx = NULL) {
  kappa <- .get_y(fit, kidx) - 0.5 * .get_N(fit, kidx)
  return(kappa)
}

#' Compute E[w] where w are the PG random variables in the augmented model
compute_omega <- function(fit, kidx = NULL) {
  omega <- pg_mean(.get_N(fit, kidx), .get_xi(fit, kidx))
  return(omega)
}

#' Compute a scaled posterior mean
#' TODO: add support for non-zero prior mean
#' @param fit a SER object
compute_nu.binser <- function(fit, idx = NULL, kidx = NULL, shift = 0) {
  prior_mean <- .get_prior_mean(fit, idx = idx)
  tau0 <- 1. / .get_var0(fit, idx = idx)

  kappa <- compute_kappa(fit, kidx)
  Zd <- compute_Zd.binser(fit, idx, shift)
  omega <- compute_omega(fit)
  tmp <- kappa - omega * Zd
  nu <- Matrix::drop(tmp %*% fit$data$X) + (prior_mean * tau0)
  return(nu)
}

#' Compute partial posterior variance
#' use this immediately after updating xi
compute_tau <- function(fit) {
  omega <- compute_omega(fit)
  tau <- Matrix::drop(omega %*% fit$data$X2)
  return(tau)
}


#' update for variational parameter parameter xi q(w) = PG(N, xi)
update_xi.binser <- function(fit, shift = 0) {
  Xb2 <- compute_Xb2.binser(fit, shift = shift)
  xi <- sqrt(abs(Xb2))
  return(xi)
}

#' update intercept and fixed effect covariates
update_delta.binser <- function(fit, idx = NULL, kidx = NULL, shift = 0) {
  Z <- fit$data$Z
  omega <- compute_omega(fit)
  Xb <- fit$data$X %*% coef.binser(fit, idx) + shift
  kappa <- compute_kappa(fit, kidx)
  delta <- Matrix::drop(Matrix::solve((omega * t(Z)) %*% Z, t(Z) %*% (kappa - omega * Xb)))
  return(delta)
}

#' update SER variational parameters
#' @param nu
#' @param tau
#' @param pi prior probability vector (length p)
.update_b_ser <- function(nu, tau, pi) {
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
update_b.binser <- function(fit, idx = NULL, kidx = NULL, shift = 0) {
  nu <- Matrix::drop(compute_nu.binser(fit, idx, kidx, shift))
  tau0 <- 1 / .get_var0(fit, idx)
  tau <- tau0 + fit$params$tau
  post <- .update_b_ser(nu, tau, .get_pi(fit, idx))
  return(post)
}


#' Variational M step for updating prior_variance
#' Optimize ELBO wrt to prior effect standard deviation
update_prior_variance.binser <- function(fit, idx = NULL) {
  b2 <- .get_alpha(fit, idx) * (.get_mu(fit, idx)^2 + .get_var(fit, idx))
  b <- coef.binser(fit, idx)
  prior_mean <- .get_prior_mean(fit, idx)
  prior_variance <- sum(b2 - 2 * b * prior_mean) + prior_mean^2
  return(prior_variance)
}


#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
#' @returns a vector of length N with the lower bound for each data point
jj_bound.binser <- function(fit, shift = 0, kidx = NULL) {
  xi <- .get_xi(fit, kidx = kidx)
  Xb <- compute_Xb.binser(fit, shift = shift)
  kappa <- compute_kappa(fit, kidx = kidx)
  n <- .get_N(fit)

  # Xb2 <- compute_Xb2.binser(fit)
  # omega <- compute_omega(fit)

  bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi) #+ 0.5 * omega * (xi^2 - Xb2)
  return(bound)
}

explicit_elbo.binser <- function(fit) {
  Xb <- compute_Xb.binser(fit)
  kappa <- compute_kappa(fit)
  xi <- .get_xi(fit)
  n <- .get_N(fit)

  Xb2 <- compute_Xb2.binser(fit)
  omega <- compute_omega(fit)

  ll <- -n * log(2) + kappa * Xb - 0.5 * Xb2 * omega
  kl <- pg_kl(n, xi)

  return(ll - kl)
}

#' More explicit, only works for Bernoulli observations
jj_bound.logistic <- function(fit, shift = 0, kidx = NULL) {
  xi <- .get_xi(fit, kidx)
  Xb <- compute_Xb.binser(fit, shift = shift)
  kappa <- compute_kappa(fit, kidx)

  Xb2 <- compute_Xb2.binser(fit, shift = shift)
  omega <- compute_omega(fit)
  bound <- kappa * Xb - 0.5 * xi - 0.25 / xi * tanh(xi / 2) * (Xb2 - xi^2) - log(1 + exp(-xi))
  return(bound)
}


compute_elbo.binser <- function(fit, shift = 0) {
  jj <- sum(jj_bound.binser(fit, shift = shift)) # E[p(y | w, b)] - KL[q(b) || p(b)]
  kl <- compute_kl.binser(fit) # KL[q(b) || p(b)]
  return(jj - kl)
}


#' initialize SER
#' @param data a list containing X, Z, y
#' @param prior_mean prior effect mean
#' @param simga0 prior effect standard deviation
init.binser <- function(data, prior_mean = 0, prior_variance = 1, prior_weights = NULL) {
  n <- nrow(data$X)
  p <- ncol(data$X) # number of covariates to select
  p2 <- ncol(data$Z) # number of fixed effects + intercept

  # initialize variational parameters
  params <- list(
    alpha = rep(1, p) / p,
    mu = rep(0, p),
    var = rep(1, p),
    delta = rep(0, p2),
    xi = rep(1e-3, n),
    tau = 0 # initialized at end of this function
  )

  # initialize hyper-parameters
  if (is.null(prior_weights)) {
    prior_weights <- rep(1, p) / p
  }
  hypers <- list(
    pi = prior_weights, # categorical probability of nonzero effect
    prior_mean = prior_mean, # prior mean of effect
    prior_variance = prior_variance # prior variance of effect, TODO: include as init arg
  )

  data$X2 <- data$X^2

  # TODO: check that data has y, N, X, Z
  fit.init <- list(
    data = data,
    params = params,
    hypers = hypers,
    elbo = c(-Inf)
  )
  fit.init$params$tau <- compute_tau(fit.init)
  return(fit.init)
}

iter.binser <- function(fit, shift = 0, fit_intercept = TRUE, fit_prior_variance = TRUE) {
  # update b
  post <- update_b.binser(fit, shift = shift)
  fit$params$mu <- post$mu
  fit$params$var <- post$var
  fit$params$alpha <- post$alpha

  # update intercept/fixed effect covariates
  if (fit_intercept) {
    fit$params$delta <- update_delta.binser(fit, shift = shift)
  }
  if (fit_prior_variance) {
    fit$hypers$prior_variance <- update_prior_variance.binser(fit)
  }

  # update xi
  fit$params$xi <- update_xi.binser(fit, shift = shift)
  fit$params$tau <- compute_tau(fit)

  return(fit)
}

#' Fit the binomial single effect regression
fit.binser <- function(data = NULL, fit = NULL,
                       maxiter = 10,
                       tol = 1e-3,
                       fit_intercept = TRUE,
                       fit_prior_variance = TRUE,
                       shift = NULL) {
  # no shift
  if (is.null(shift)) {
    shift <- 0
  }

  # initialize model if no provided
  if (is.null(fit)) {
    fit <- init.binser(data)
  }

  # iteration loop
  for (i in 1:maxiter) {
    fit <- iter.binser(fit)

    # compute elbo -- TODO: do we want to record every time?
    fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
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
  # Initialize model
  data <- binsusie_prep_data(X, y, N, Z, scale = scale, center = center)
  if (is.null(s_init)) {
    fit <- init.binser(data, prior_mean = prior_mean, prior_variance = prior_variance, prior_weights = prior_weights)
  } else {
    fit <- s_init
  }

  # model fitting
  for (i in 1:max_iter) {
    fit <- iter.binser(
      fit,
      shift = o,
      fit_intercept = intercept,
      fit_prior_variance = estimate_prior_variance
    )
    fit$elbo <- c(fit$elbo, compute_elbo.binser(fit, shift = o))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
}
