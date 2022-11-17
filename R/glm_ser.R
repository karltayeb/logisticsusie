# Adapted from: https://andrewg3311.github.io/susieR_logistic_wflow/susie_logistic_demonstration.html#adjustments_to_susie_code

# compute SER using GLM
#' @export
fit_glm_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, prior_weights = NULL, family = "binomial") {
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
  return(list(alpha = alpha, mu = post_mean, intercept = intercept, var = post_var, lbf = lbf, lbf_model = lbf_model, prior_variance = prior_variance, loglik = loglik))
}
