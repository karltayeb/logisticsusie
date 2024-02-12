#' Take a function for fitting MLE of univariate regression and make an SER function
#'
#' this is a utlity function for generating SER fit function compatible with `ibss_from_ser`
#' if you have a function for fitting a univariate regression this will create a function for fitting an SER
#' it relied on an asymptotic approximation-- so if you have betahat, standard error, and likelihood ratio you can use this
#'
#' @param uni_fun a function for fitting univariate regression
#'  must take arguments: x, y, o, estimate_intercept
#'    additional arguments passed via ...
#'  must return: a list with items lbf = log bayes factor, mu = posterior mean, var = posterior variance, intercept=intercept
#'    note that posterior variance doesn't figure in computations so you can set it to a dummy value (e.g. 0)
#' @export
ser_from_univariate <- function(uni_fun) {
  # a function for fitting SER
  ser_fun <- function(X, y, o = NULL, prior_variance = 1, estimate_prior_variance=T, estimate_intercept = T, prior_weights = NULL, ...) {
    # set=up
    p <- dim(X)[2]

    # 0 offset if not specified
    if (is.null(o)) {
      o <- rep(0, length(y))
    }

    # use uvb intercepts as fixed intercept
    fits <- dplyr::bind_rows(purrr::map(1:p, ~ uni_fun(X[, .x], y, o, prior_variance = prior_variance, estimate_intercept = estimate_intercept, ...)))

    # compute summaries
    alpha <- exp(fits$lbf - matrixStats::logSumExp(fits$lbf))
    lbf_model <- sum(alpha * fits$lbf) - categorical_kl(alpha, rep(1 / p, p))

    # return standard ouput: mu var alpha intercept, lbf
    res <- list(
      mu = fits$mu,
      var = fits$var,
      alpha = alpha,
      intercept = fits$intercept,
      lbf = fits$lbf,
      lbf_model = lbf_model,
      prior_variance = prior_variance
    )
    class(res) <- "ser"
    return(res)
  }
  return(ser_fun)
}


#' Compute log of the asymptotic Bayes factor (ABF) for the SER vs null model
#'
#' @param betahat vector of MLEs for the effect size
#' @param s2 vector of standard errors
#' @param prior_variance variance of the effect size b ~ N(0, prior_variance)
#' @returns log ABF for each variable
compute_log_abf <- function(betahat, s2, prior_variance){
  log_abf <- dnorm(betahat, 0, sd = sqrt(s2 + prior_variance), log = T) -
    dnorm(betahat, 0, sd = sqrt(s2), log = T)
  return(log_abf)
}

#' Compute log of the asymptotic Bayes factor (ABF) for the SER vs null model
#'
#' @param betahat vector of MLEs for the effect size
#' @param s2 vector of standard errors
#' @param prior_variance variance of the effect size b ~ N(0, prior_variance)
#' @param pi prior over non-zero effects
#' @returns log of the ABF for the SER log(\sum_j \pi_j ABF_j)
compute_log_abf_ser <- function(betahat, s2, prior_variance, pi){
  log_abf <- compute_log_abf(betahat, s2, prior_variance)
  return(logSumExp(log_abf + log(pi)))
}


#' Compute log of the CORRECTED asymptotic Bayes factor (ABF)
#'
#' @param betahat vector of MLEs for the effect size
#' @param s2 vector of standard errors
#' @param prior_variance variance of the effect size b ~ N(0, prior_variance)
#' @returns log of the LABF for each variable
compute_log_labf <- function(betahat, s2, lr, prior_variance){
  log_abf <- compute_log_abf(betahat, s2, prior_variance)
  alr_mle <- dnorm(betahat, betahat, sqrt(s2), log=T) -
    dnorm(betahat, 0, sqrt(s2), log=T)
  log_labf <- log_abf - alr_mle + lr
  return(log_labf)
}

#' Compute log of the CORRECTED asymptotic Bayes factor (ABF) for the SER vs null model
#'
#' @param betahat vector of MLEs for the effect size
#' @param s2 vector of standard errors
#' @param prior_variance variance of the effect size b ~ N(0, prior_variance)
#' @param pi prior over non-zero effects
#' @returns log of the ABF for the SER log(\sum_j \pi_j ABF_j)
compute_log_labf_ser <- function(betahat, s2, lr, prior_variance, pi){
  log_labf <- compute_log_labf(betahat, s2, lr, prior_variance)
  return(logSumExp(log_labf + log(pi)))
}

#' Optimize prior variance
#'
#' Optimize the asymptotic BF for the SER w.r.t the prior variance
#' @param betahat vector of MLEs for the effect size
#' @param s2 vector of standard errors
#' @param pi prior over non-zero effects
#' @returns the optimimum prior variance
optimize_prior_variance <- function(betahat, s2, lr, pi, laplace=T, min_prior_variance=0){
  if(laplace){
    f_opt <- function(prior_variance){
      compute_log_labf_ser(betahat, s2, lr, prior_variance, pi)
    }
  } else{
    f_opt <- function(prior_variance){
      compute_log_abf_ser(betahat, s2, prior_variance, pi)
    }
  }

  # TODO: reasonable upper bound?
  upper <- max(4 * betahat^2)
  opt <- optimise(f_opt, c(0, upper), maximum = T)$maximum

  fopt <- f_opt(opt)
  f0 <- f_opt(0)

  # set prior variance to 0 if optimized results is not better than null
  if(f0 >= fopt){
    opt <- 0
  }

  # set prior variance to 0 if it is not bigger than min_prior_variance
  if(opt <= min_prior_variance){
    opt <- 0
  }
  return(opt)
}

asymptotic_ser <- function(betahat, shat2, intercept, lr, prior_weights, prior_variance = 1, estimate_prior_variance=T, laplace=T, min_prior_variance=0){
  p <- length(betahat)

  # estimate prior variance
  if(estimate_prior_variance){
    prior_variance <- optimize_prior_variance(
      betahat, shat2, lr, prior_weights, laplace, min_prior_variance)
  }

  # compute ABF
  log_abf <- dnorm(betahat, 0, sqrt(prior_variance + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)

  # compute liklihood ratio under asymptotic approximation
  alr_mle <- dnorm(betahat, betahat, sqrt(shat2), log=T) -
    dnorm(betahat, 0, sqrt(shat2), log=T)

  # make corrected ABF
  if(laplace){
    lbf <- log_abf - alr_mle + lr
  } else{
    lbf <- log_abf
  }

  # deal with perfectly enriched/depleted
  lbf[is.infinite(shat2)] <- 0

  maxlbf <- max(lbf)
  w <- exp(lbf - maxlbf)
  w_weighted <- w * prior_weights
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted/weighted_sum_w

  # posterior
  post_var <- 1/((1/shat2) + (1/prior_variance))
  if(prior_variance != 0){
    post_mean <- (1/shat2) * post_var * betahat
  } else {
    # if prior variance is 0, posterior mean is 0
    post_mean <- rep(0, p)
  }
  post_mean2 <- post_var + post_mean^2
  lbf_model <- maxlbf + log(weighted_sum_w)

  res <- list(
    alpha = alpha,
    mu = post_mean,
    var = post_var,
    lbf = lbf,
    log_abf = log_abf,
    lr = lr,
    alr_mle = alr_mle,
    lbf_model = lbf_model,
    prior_variance = prior_variance,
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept
  )
  return(res)
}

#' Take a function for fitting univariate regression and make an SER function
#'
#' this is a utlity function for generating SER fit function compatible with `ibss_from_ser`
#' if you have a function for fitting a univariate regression this will create a function for fitting an SER
#'
#' @param uni_fun a function for fitting univariate regression
#'  must take arguments: x, y, o, prior_variance, estimate_intercept
#'    additional arguments passed via ...
#'    hint: it may be useful to pass null model likelihood in ...
#'  must return: a list with items betahat = effect estimate, shat2 = square of standard errror, lr = likelihood ratio, intercept=intercept
#' @export
asymptotic_ser_from_univariate <- function(uni_fun) {
  # a function for fitting SER
  ser_fun <- function(X, y, o = NULL, prior_variance = 1, estimate_prior_variance=T, estimate_intercept = T, prior_weights = NULL, laplace = T, min_prior_variance=0., ...) {
    # set=up
    p <- dim(X)[2]

    # 0 offset if not specified
    if (is.null(o)) {
      o <- rep(0, length(y))
    }

    # default pi
    if (is.null(prior_weights)) {
      prior_weights <- rep(1/p, p)
    }

    # fit mles
    uni <- purrr::map(1:p, ~ uni_fun(X[, .x], y, o, estimate_intercept = estimate_intercept, ...))

    # extract relevant values
    lr <- purrr::map_dbl(uni, ~purrr::pluck(.x, 'lr'))
    betahat <- purrr::map_dbl(uni, ~purrr::pluck(.x, 'betahat'))
    shat2 <- purrr::map_dbl(uni, ~purrr::pluck(.x, 'shat2'))
    intercept <- purrr::map_dbl(uni, ~purrr::pluck(.x, 'intercept'))

    res <- asymptotic_ser(betahat, shat2, intercept, lr, prior_weights, estimate_prior_variance, laplace, min_prior_variance)
    return(res)
  }
  return(ser_fun)
}

