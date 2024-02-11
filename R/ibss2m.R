# Heuristic Binomial SuSiE, iteratively fit SERs
# but SERs don't need to share variational parameters

# Computations ----

###
# Updates ----
###

#' update alpha, mu, and var
update_model.ibss2m <- function(fit, data,
                               fit_intercept = TRUE,
                               fit_prior_variance = TRUE,
                               fit_alpha = TRUE,
                               fit_xi = TRUE,
                               update_idx = NULL,
                               track_elbo=TRUE) {
  if (is.null(update_idx)) {
    update_idx <- seq(fit$L)
  }

  # on first iteration we don't have re stored... compute from scratch
  if(is.null(fit$re)){
    re <- list(mu = 0, mu2 = 0)
    add_only <- T
  } else{
    re <- fit$re
    add_only <- F
  }

  for (l in update_idx) {
    # remove current effect estimate
    if(!add_only){
      re = sub_re(re$mu,
                  re$mu2,
                  compute_Xb(fit$sers[[l]], data),
                  compute_Xb2(fit$sers[[l]], data))
    }

    data$shift <- re$mu
    data$shift_var <- re$mu2 - re$mu^2

    # update SER
    fit$sers[[l]] <- update_model(fit$sers[[l]], data,
                                  fit_intercept = fit_intercept,
                                  fit_prior_variance = fit_prior_variance)


    # add current effect estimate
    data$shift <- 0
    data$shift_var <- 0
    re = add_re(re$mu,
                re$mu2,
                compute_Xb(fit$sers[[l]], data),
                compute_Xb2(fit$sers[[l]], data))
  }

  # NOTE: there is no elbo for this
  # but the sum of the elbo of SERs should converge when the model converes
  elbo <- sum(purrr::map_dbl(fit$sers, ~tail(.x$elbo, 1)))
  if(track_elbo){
    fit$elbo <- c(fit$elbo, elbo)
  }
  fit$re <- re
  return(fit)
}

# Initialization -------
initialize_ibss2m <- function(L, n, p, p2, mu0=0, var0=1, pi = rep(1/p, p)){
  if(length(mu0) < L){
    mu0 <- rep(mu0, L)
  }
  if(length(var0) < L){
    var0 <- rep(var0, L)
  }
  sers <- purrr::map(1:L, ~initialize_uvbser(n, p, p2, mu0[.x], var0[.x], pi))

  # set prior variance to 0 so that compute_psi2 = 0-- better initial iteration
  # for(l in 1:L){
  #   sers[[l]]$var <- rep(0, p)
  # }

  fit <- list(sers=sers, L=L, elbo=-Inf)
  class(fit) <- c('ibss2m', 'binsusie')
  return(fit)
}

data_initialize_ibss2m <- function(data, L, mu0=0, var0=1){
  n <- nrow(data$X)
  p <- ncol(data$X)
  p2 <- ncol(data$Z)
  fit <- initialize_ibss2m(L, n, p, p2, mu0, var0)
  return(fit)
}

#' IBSS2M
#'
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
#'
#' @export
ibss2m <- function(X,
                     y,
                     N = rep(1, length(y)), # number of trials for each
                     Z = NULL,
                     L = min(10, ncol(X)),
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
                     tol = 0.001) {

  data <- binsusie_prep_data(X, y, N, Z, center=center, scale=scale)
  fit <- data_initialize_ibss2m(data, L, mu0 = prior_mean, var0 = prior_variance)
  fit <- fit_model(fit, data, fit_prior_variance=estimate_prior_variance, max_iter = max_iter)
  return(fit)
}
