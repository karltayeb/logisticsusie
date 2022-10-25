
#' update alpha, mu, and var
update_b_hard_select.binsusie <- function(fit, fit_intercept = TRUE, fit_prior_variance = TRUE, update_idx = NULL) {
  if (is.null(update_idx)) {
    update_idx <- seq(fit$hypers$L)
  }
  shift <- compute_Xb.binsusie(fit)
  for (l in update_idx) {
    # remove current effect estimate
    shift <- Matrix::drop(shift - compute_Xb.binser(fit, idx = l))

    for (i in 1:10) {
      # update SER
      post_l <- update_b.binser(fit, idx = l, shift = shift)
      fit$params$mu[l, ] <- post_l$mu
      fit$params$var[l, ] <- post_l$var

      ### THIS IS THE INTERVENTION -- WILL HURT ELBO
      alpha_hard <- rep(0, length(post_l$alpha))
      alpha_hard[which.max(post_l$alpha)] <- 1
      fit$params$alpha[l, ] <- alpha_hard

      # update intercept/fixed effect covariates
      if (fit_intercept) {
        fit$params$delta[l, ] <- update_delta.binser(fit, idx = l, shift = shift)
      }

      # update prior_variance
      if (fit_prior_variance) {
        fit$hypers$prior_variance[l] <- update_prior_variance.binser(fit, idx = l)
      }

      # update xi and tau
      fit$params$xi <- update_xi.binsusie(fit)
      fit$params$tau <- compute_tau(fit)
    }

    # add current effect estimate
    shift <- shift + compute_Xb.binser(fit, idx = l)
  }
  return(fit)
}


#' initialize SER
.forward_init <- function(data, L = 5, init_iter = 10, prior_mean = 0, prior_variance = 1, prior_weights = NULL, kidx = NULL) {
  fit <- init.binsusie(data, L, prior_mean, prior_variance, prior_weights, kidx)

  # update each single effect
  # but intervene by setting alpha to a hard selection
  for (i in 1:init_iter) {
    fit <- update_b_hard_select.binsusie(fit, fit_intercept = T, fit_prior_variance = F)
  }
  return(fit)
}


binsusie_forward_init <- function(X,
                                  y,
                                  N = rep(1, length(y)), # number of trials for each
                                  Z = NULL,
                                  L = min(10, ncol(X)),
                                  prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
                                  prior_variance = 1.0, # TODO: make scaled prior variance?
                                  prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect...
                                  init_iter = 5) {
  ##################
  # checking / setup

  # Make data list
  data <- list(X = X, Z = Z, y = y, N = N)

  # If N was passed as a scalar, convert to vector of length n
  if (length(data$N) == 1) {
    data$N <- rep(N, length(y))
  }
  .check_data_binsusie(data) # throws errors if something is wrong

  fit <- .forward_init(
    data,
    L = L,
    init_iter = init_iter,
    prior_mean = prior_mean,
    prior_variance = prior_variance,
    prior_weights = prior_weights
  )
  return(fit)
}


binsusie_init <- function(X,
                          y,
                          N = rep(1, length(y)), # number of trials for each
                          Z = NULL,
                          L = min(10, ncol(X)),
                          prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
                          prior_variance = 1.0, # TODO: make scaled prior variance?
                          prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect...
                          init_iter = 5) {
  ##################
  # checking / setup

  # Make data list
  data <- list(X = X, Z = Z, y = y, N = N)

  # If N was passed as a scalar, convert to vector of length n
  if (length(data$N) == 1) {
    data$N <- rep(N, length(y))
  }
  .check_data_binsusie(data) # throws errors if something is wrong

  fit <- init.binsusie(
    data,
    L = L,
    prior_mean = prior_mean,
    prior_variance = prior_variance,
    prior_weights = prior_weights
  )
  return(fit)
}
