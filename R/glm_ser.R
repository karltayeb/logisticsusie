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
map_univariate_regression <- function(X, y, off = NULL, estimate_intercept=T, family='binomial', augment=F){
  if(augment==T){
    warning('`augment == T` not implemeneted')
  }
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


#' Prepare data for fast glm
#'
#' Helper function for preparing the design matrix, response vector, and offset
#' @param x the variable of interest
#' @param y the response vector
#' @param o fixed offset, length(o) == length(y)
#' @param intercept boolean to include column of 1s in design matrix
#' @param augment boolean to augment data to ensure existence of MLE
#' @returns a list containing design `X`, reponse `y`, and offset `o`
prep_data_for_fastglm <- function(x, y, o, intercept=T, augment=F){
  if(augment){
    n <- length(y) + 4
    x <- c(x, c(1, 1, 0, 0))
    y <- c(y, c(1, 0, 1, 0))
    o <- c(o, rep(0, 4))
    weights <- c(rep(1, n-4), rep(1e-6, 4))
  } else{
    n <- length(y)
    weights <- rep(1, n)
  }

  if(intercept){
    ones <- rep(1, n)
    X <- cbind(ones, x)
  } else{
    X <- as.matrix(x, nrow=n)
  }
  regression_data <- list(X=X, y=y, o=o, weights=weights)
  return(regression_data)
}

#' Fit fast GLM
#'
#' Fit fast GLM, and extract necessary information for SER
#' @param x the variable of interest
#' @param y the response vector
#' @param o fixed offset, length(o) == length(y)
#' @param family the family of the glm, e.g. "binomial", "poisson"
#' @param augment boolean to augment data to ensure existence of MLE
#' @param fit0
fit_fast_glm <- function(x, y, o, family='binomial', augment=F, fit0=NULL){
  # 1. data setup
  regression_data <- prep_data_for_fastglm(x, y, o, intercept=T, augment=augment)

  # 2. fit model
  fit <- fastglm::fastglm(y = regression_data$y,
                          x = regression_data$X,
                          offset = regression_data$o,
                          family = family,
                          weights = regression_data$wei)

  # 3. fit null model if not provided
  if(is.null(fit0)){
    n <- length(regression_data$y)
    ones <- rep(1, n)
    fit0 <- fastglm::fastglm(y = regression_data$y,
                             x = matrix(ones, nrow = n),
                             offset = regression_data$o,
                             family = family)
  }

  # 3. processing
  coef <- unname(fit$coefficients)
  intercept <- coef[1]
  betahat <- coef[2]
  shat2 <-  fit$se[2]^2
  lr  <- 0.5 * (fit0$deviance - fit$deviance)
  res <- list(
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept,
    lr = lr
  )
  return(res)
}

#' Map univariate regression fit with fastglm
#'
#' Fit univariate glm for each column of X
#' @param X n x p design matrix
#' @param y n vector response
#' @param o optional offset, vector of length n
#' @param estimate_intercept boolean to fit intercept (not implemented)
#' @param family family for glm, e.g. binomial
#' @returns list containing intercept and effect estimates,
#'    standard errors, and likelihood ratios
map_fastglm <- function(X, y, off = NULL, estimate_intercept=T, family='binomial', augment=T){
  if(!estimate_intercept){
    warning("estimate_intercept = FALSE not implemented")
  }
  p <- ncol(X)
  if (is.null(off)) {
    off <- rep(0, length(y))
  }

  # 1. fit null model once
  null_data <- prep_data_for_fastglm(x = rep(1, length(y)), y = y, o = off,
                                     intercept = F, augment = augment)
  fit0 <- fastglm::fastglm(y = null_data$y,
                           x = null_data$X,
                           offset = null_data$o,
                           family = family,
                           weights = null_data$weights)

  # 2. fit each variable
  res <- purrr::map_dfr(1:p, ~fit_fast_glm(X[, .x], y, off, family=family, augment = augment, fit0 = fit0))

  # 3. process results
  res <- list(
    p = p,
    betahat = res$betahat,
    shat2 = res$shat2,
    intercept = res$intercept,
    lr = res$lr
  )
  return(res)
}

#' compute SER using GLM, use Laplace approximation to BF rather than ABF
#' @export
fit_glm_ser2 <- function(X, y, o = NULL,
                         prior_variance = 1, intercept = T,
                         prior_weights = NULL, family = "binomial",
                         laplace = T, estimate_prior_variance=T,
                         augment = T,
                         glm_mapper = map_fastglm) {

  # univariate regressions
  uni <- glm_mapper(X, y, o,
                    estimate_intercept=intercept,
                    family=family,
                    augment=augment)
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
    log_abf = log_abf,
    lr = lr,
    alr_mle = alr_mle,
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
