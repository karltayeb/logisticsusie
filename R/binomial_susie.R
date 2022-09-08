####
# Binomial SuSiE
###


###
# Computations
##


#' Get Binomial SER coefficients
#' Extract regression coefficients from binomial SER fit
#' @param fit Binomial SER object
#' @return Return E[\beta]
coef.binsusie <- function(fit, idx=NULL){
  b <- colSums(.get_alpha(fit, idx) * .get_mu(fit, idx))
  return(b)
}

#' Expected linear prediction
compute_Xb.binsusie <- function(fit){
  Xb <- drop(fit$data$X %*% colSums(.get_alpha(fit) * .get_mu(fit)))
  Zd <- drop(fit$data$Z %*% colSums(.get_delta(fit)))
  return(Xb + Zd)
}


#' Second moment of linear prediction
compute_Xb2.binsusie <- function(fit){
  B <- .get_alpha(fit) * .get_mu(fit)
  XB <- fit$data$X %*% t(B)
  Xb <- rowSums(XB)
  Zd <- drop(fit$data$Z %*% colSums(.get_delta(fit)))

  B2 <- .get_alpha(fit) * (.get_mu(fit)^2 + .get_var(fit))
  b2 <- colSums(B2)

  Xb2 <- fit$data$X^2 %*% b2 + Xb**2 - rowSums(XB^2)
  Xb2 <- drop(Xb2 + 2*Xb*Zd + Zd^2)
  return(Xb2)
}

#' Compute KL(q(\beta) || p(\beta)) for Sum of SERs
#' Note this does not include KL(q(\omega) || p(\omega))
#' since that gets computes as part of the JJ bound
compute_kl.binsusie <- function(fit){
  kl <- sum(purrr::map_dbl(seq(fit$hypers$L), ~ .compute_ser_kl(
    .get_alpha(fit, .x),
    .get_pi(fit, .x),
    .get_mu(fit, .x),
    .get_var(fit, .x),
    .get_mu0(fit, .x),
    .get_var0(fit, .x)
  )))
  return(kl)
}


#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
jj_bound.binsusie <- function(fit){
  xi <- .get_xi(fit)
  Xb <- compute_Xb.binsusie(fit)
  kappa <- compute_kappa(fit)
  n <- .get_N(fit)

  Xb2 <- compute_Xb2.binsusie(fit)
  omega <- compute_omega(fit)

  bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi) #+ 0.5 * omega * (xi^2 - Xb2)
  return(bound)
}


compute_elbo.binsusie <- function(fit){
  jj <- sum(jj_bound.binsusie(fit))
  kl <- compute_kl.binsusie(fit)
  return(jj - kl)
}

###
# Updates
###


#' update alpha, mu, and var
update_b.binsusie <- function(fit, fit_intercept=TRUE, fit_prior_variance=TRUE){
  shift <- compute_Xb.binsusie(fit)
  for(l in seq(fit$hypers$L)){
    # remove current effect estimate
    shift <- shift - compute_Xb.binser(fit, idx=l)

    # update SER
    post_l <- update_b.binser(fit, idx=l, shift=shift)
    fit$params$mu[l,] <- post_l$mu
    fit$params$var[l,] <- post_l$var
    fit$params$alpha[l,] <- post_l$alpha

    # update intercept/fixed effect covariates
    if(fit_intercept){
      fit$params$delta[l,] <- update_delta.binser(fit, idx=l, shift=shift)
    }

    # update sigma0
    if(fit_prior_variance){
      fit$hypers$sigma0[l] <- update_sigma0.binser(fit, idx=l)
    }

    # add current effect estimate
    shift <- shift + compute_Xb.binser(fit, idx=l)
  }
  return(fit)
}


#' update for variational parameter parameter xi q(w) = PG(N, xi)
update_xi.binsusie <- function(fit){
  Xb2 <- compute_Xb2.binsusie(fit)
  xi <- sqrt(abs(Xb2))
  return(xi)
}

###
# Fitting
###

.init.binsusie.params <- function(n, p, p2, L){
  params <- list(
    alpha = matrix(rep(1, p*L)/p, nrow=L),
    mu = matrix(rep(0, L*p), nrow=L),
    var = matrix(rep(1, L*p), nrow=L),
    delta = matrix(rep(0, p2*L, nrow=L)),
    xi = rep(1e-3, n),
    tau = 0 # initialized at end of this function
  )
  return(params)
}

.init.binsusie.hypers <- function(n, p, L, prior_mean, prior_variance, prior_weights){

  stopifnot("`prior_mean` must be scalar of length L`"=length(prior_mean) %in% c(1, L))
  if(length(prior_mean) == 1){
    prior_mean <- rep(prior_mean, L)
  }
  stopifnot("`prior_variance` must be scalar of length L`"=length(prior_variance) %in% c(1, L))
  if(length(prior_variance) == 1){
    prior_variance <- rep(prior_variance, L)
  }

  if(is.null(prior_weights)){
    prior_weights <- matrix(rep(1/p, p*L), nrow=L)
  } else if(is.matrix(prior_weights)){
    stopifnot("`prior weights must be a L x p matrix or p vector" = all(dim(prior_weights) == c(L, p)))
    stopifnot("`prior weights must sum to one" = all(rowSums(prior_weights)==1))
  } else{
    stopifnot("`prior weights must be a L x p matrix or p vector" = length(prior_weights) == p)
    stopifnot("`prior weights must sum to one" = sum(prior_weights) == 1)
    prior_weights <- matrix(rep(prior_weights, L), nrow = L, byrow = T)
  }

  hypers <- list(
    L = L,
    pi = prior_weights,  # categorical probability of nonzero effect
    mu0 = prior_mean,  # prior mean of effect
    sigma0 = prior_variance  # prior variance of effect, TODO: include as init arg,
  )
  return(hypers)
}

#' initialize SER
init.binsusie <- function(
    data, L=5, prior_mean=0, prior_variance = 1, prior_weights=NULL, kidx=NULL){

  n <- nrow(data$X)
  p <- ncol(data$X)

  # set Z if not provided
  if(is.null(data$Z)){
    data$Z <- matrix(rep(1, n), nrow = n)
  }
  p2 <- ncol(data$Z)

  params <- .init.binsusie.params(n, p, p2, L)
  hypers <- .init.binsusie.hypers(n, p, L, prior_mean, prior_variance, prior_weights)

  # TODO: check that data has y, N, X, Z
  fit.init <- list(
    data = data,
    params = params,
    hypers = hypers,
    elbo = c(-Inf)
  )
  fit.init$kidx <- kidx
  fit.init$params$tau <- compute_tau(fit.init)
  return(fit.init)
}


iter.binsusie <- function(fit, fit_intercept=TRUE, fit_prior_variance=TRUE){
  # update b
  fit <- update_b.binsusie(fit, fit_intercept, fit_prior_variance)
  # update xi
  fit$params$xi <- update_xi.binsusie(fit)
  fit$params$tau <- compute_tau(fit)
  return(fit)
}

#' Fit the binomial single effect regression
fit.binsusie <- function(
    data,
    maxiter=10,
    tol=1e-3,
    fit_intercept=TRUE,
    fit_prior_variance=TRUE, kidx=NULL){

  fit <- init.binsusie(data, kidx=kidx)

  for(i in 1:maxiter){
    fit <- iter.binsusie(fit, fit_intercept, fit_prior_variance)
    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.binsusie(fit))
    if (.converged(fit, tol)){
      break
    }
  }
  return(fit)
}

#' @param X a n x p matrix of covariates
#' @param y an n vector of integer counts, bernoulli/binomial observations
#' @param N the number of binomial trials, defaults to 1, may be a scalar or vector of length n
#' @param Z fixed effect covaraites (including intercept). If null just a n x 1 matrix of ones
#' @param prior_mean the prior mean of each non-zero element of b. Either a scalar or vector of length L.
#' @param prior_variance the prior variance of each non-zero element of b. Either a scalar or vector of length L. If `estimate_prior_variance=TRUE` the value provides an initial estimate of prior variances
#' @param prior_weights prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix
#' @param intercept
#' @param estimate_prior_variance
#' @param s_init a logistic susie object to initialize with, NOTE if non-null, we ignore `prior_mean`, `prior_variance`, and `prior_weights`
#' @export
binsusie <- function(
  X,
  y,
  N = 1,  # number of trials for each
  Z=NULL,
  L = min(10, ncol(X)),
  prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
  prior_variance = 1.0,  # TODO: make scaled prior variance?
  prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect... = hypers$pi
  #null_weight = NULL,
  #standardize = TRUE,
  intercept = TRUE,
  estimate_prior_variance = TRUE,
  #estimate_prior_method = c("optim", "EM", "simple"),  right now only EM is implimented
  #check_null_threshold = 0,
  #prior_tol = 1e-09,
  s_init = NULL, # previous fit with which to initialize NO
  coverage = 0.95,
  min_abs_corr = 0.5,
  max_iter = 100,
  tol = 0.001,
  verbose = FALSE,
  #track_fit = FALSE,
  #refine = FALSE,
  n_purity = 100
){
  ##################
  # checking / setup

  # Make data list
  data <- list(X = X, Z=Z, y=y, N=N)

  # Initialize object
  if(is.null(s_init)){
    fit <- init.binsusie(
      data,
      L=L,
      prior_mean = prior_mean,
      prior_variance = prior_variance,
      prior_weights=prior_weights
    )
  } else{
    fit <- s_init
  }

  # model fitting
  for(i in 1:max_iter){
    fit <- iter.binsusie(fit, intercept, estimate_prior_variance)
    fit$elbo <- c(fit$elbo, compute_elbo.binsusie(fit))
    if (.converged(fit, tol)){
      break
    }
  }

  # post-processing

  # 1. compute credible sets
  return(fit)
}

