# Adapted from: https://andrewg3311.github.io/susieR_logistic_wflow/susie_logistic_demonstration.html#adjustments_to_susie_code

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
  return(logSumExp(log_labf + log(pi)))
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
optimize_prior_variance <- function(betahat, s2, lr, pi, laplace=T){
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
  upper <- max(1000 * s2)
  opt <- optimise(f_opt, c(0, upper), maximum = T)$maximum
  return(opt)
}

# compute SER using GLM
#' @export
fit_glm_ser <- function(X, y, o = NULL,
                        prior_variance = 1.0, intercept = T,
                        prior_weights = NULL, family = "binomial") {
  estimate_intercept <- intercept # hack
  p <- ncol(X)
  betahat <- numeric(p)
  shat2 <- numeric(p)
  intercept <- rep(0, p)

  if (is.null(o)) {
    # fixed is components from previous SER fits
    o <- rep(0, length(y))
  }

  for (j in 1:p) {
    # logistic regression on each column of X separately
    if (estimate_intercept) {
      log.fit <- glm(y ~ X[, j] + 1 + offset(o), family = family) # fit w/ intercept
      intercept[j] <- unname(coef(log.fit)[1])
    } else {
      log.fit <- glm(y ~ X[, j] - 1 + offset(o), family = family) # fit w/out intercept
    }
    log.fit.coef <- summary(log.fit)$coefficients
    # NOTE: coerces "intercept" to be 0 or 1 to grab relevant row of glm coefficient output
    betahat[j] <- ifelse(is.na(log.fit.coef[1 + estimate_intercept, 1]), 0, log.fit.coef[1 + estimate_intercept, 1]) # beta-hat MLE (if na, just set to 0)
    shat2[j] <- ifelse(is.na(log.fit.coef[1 + estimate_intercept, 2]), Inf, log.fit.coef[1 + estimate_intercept, 2]^2) # (std errof beta-hat MLE)^2 (if na, just set to Inf)
  }

  # = wakefield ABF
  lbf <- dnorm(betahat, 0, sqrt(prior_variance + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)

  # log(bf) on each SNP
  lbf[is.infinite(shat2)] <- 0 # deal with special case of infinite shat2 (eg happens if X does not vary)


  # uniform prior on which column of X to select
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / p, p)
  }

  maxlbf <- max(lbf)
  w <- exp(lbf - maxlbf) # w is proportional to BF, but subtract max for numerical stability

  # posterior prob on each SNP
  w_weighted <- w * prior_weights
  weighted_sum_w <- sum(w_weighted)
  alpha <- w_weighted / weighted_sum_w
  post_var <- 1 / ((1 / shat2) + (1 / prior_variance)) # posterior variance
  post_mean <- (1 / shat2) * post_var * betahat # posterior mean
  post_mean2 <- post_var + post_mean^2 # posterior second moment

  # BF for single effect model
  lbf_model <- maxlbf + log(weighted_sum_w)

  # TODO: get offsets in the null model!
  null_likelihood <- dbinom(y, 1, mean(y), log = T)
  loglik <- lbf_model + null_likelihood

  res <- list(
    alpha = alpha,
    mu = post_mean,
    var = post_var,
    lbf = lbf,
    lbf_model = lbf_model,
    prior_variance = prior_variance,
    loglik = loglik,
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept
  )
  return(res)
}

#' Map univariate regression
#'
#' Fit univariate glm for each column of X
#' @param X n x p design matrix
#' @param y n vector response
#' @param o optional offset, vector of length n
#' @param estimate_intercept boolean to fit intercept
#' @param family family for glm, e.g. binomial
#' @returns list containing intercept and effect estimates,
#'    standard errors, and likelihood ratios
map_univariate_regression <- function(X, y, off = NULL, estimate_intercept=T, family='binomial'){
  p <- ncol(X)
  betahat <- numeric(p)
  shat2 <- numeric(p)
  intercept <- rep(0, p)
  lr <- rep(0, p)

  if (is.null(off)) {
    off <- rep(0, length(y))
  }

  for (j in 1:p) {
    if (estimate_intercept) {
      log.fit <- glm(y ~ X[, j] + 1 + offset(off), family = family)
      intercept[j] <- unname(coef(log.fit)[1])
    }
    else {
      log.fit <- glm(y ~ X[, j] - 1 + offset(off), family = family)
    }
    log.fit.coef <- summary(log.fit)$coefficients
    betahat[j] <- ifelse(is.na(log.fit.coef[1 + estimate_intercept,
                                            1]), 0, log.fit.coef[1 + estimate_intercept, 1])
    shat2[j] <- ifelse(is.na(log.fit.coef[1 + estimate_intercept,
                                          2]), Inf, log.fit.coef[1 + estimate_intercept, 2]^2)
    lr[j] <- 0.5 * (log.fit$null.deviance - log.fit$deviance)
  }

  res <- list(
    p = p,
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept,
    lr = lr
  )
  return(res)
}

#' compute SER using GLM, use Laplace approximation to BF rather than ABF
#' @export
fit_glm_ser2 <- function(X, y, o = NULL,
                         prior_variance = 1, intercept = T,
                         prior_weights = NULL, family = "binomial",
                         laplace=T, estimate_prior_variance=T) {
  estimate_intercept <- intercept # hack

  # univariate regressions
  uni <- map_univariate_regression(X, y, o, intercept, family)
  lr <- uni$lr
  betahat <- uni$betahat
  shat2 <- uni$shat2
  intercept <- uni$intercept
  p <- ncol(X)

  # default pi
  if (is.null(prior_weights)) {
    prior_weights <- rep(1/p, p)
  }

  # estimate prior variance
  if(estimate_prior_variance){
    prior_variance <- optimize_prior_variance(
      betahat, shat2, lr, prior_weights, laplace)
  }

  # compute ABF
  lbf <- dnorm(betahat, 0, sqrt(prior_variance + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)

  # make corrected ABF
  if(laplace){
    # compute liklihood ratio under asymptotic approximation
    alr_mle <- dnorm(betahat, betahat, sqrt(shat2), log=T) -
      dnorm(betahat, 0, sqrt(shat2), log=T)

    lbf <- lbf - alr_mle + lr
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
  post_mean <- (1/shat2) * post_var * betahat
  post_mean2 <- post_var + post_mean^2
  lbf_model <- maxlbf + log(weighted_sum_w)
  null_loglik <- sum(dbinom(y, 1, mean(y), log = T))
  loglik <- lbf_model + null_loglik

  res <- list(
    alpha = alpha,
    mu = post_mean,
    var = post_var,
    lbf = lbf,
    lr = lr,
    null_loglik = null_loglik,
    lbf_model = lbf_model,
    prior_variance = prior_variance,
    loglik = loglik,
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept
  )
  return(res)
}
