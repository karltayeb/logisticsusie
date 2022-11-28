# Adapted from: https://andrewg3311.github.io/susieR_logistic_wflow/susie_logistic_demonstration.html#adjustments_to_susie_code

# fit GLM for each column of X
map_glm <- function(X, y, o = NULL, estimate_intercept = T, family = "binomial") {
  # setup
  p <- ncol(X)
  betahat <- numeric(p)
  shat2 <- numeric(p)
  intercept <- rep(0, p)

  # 0 offset if not specified
  if (is.null(o)) {
    o <- rep(0, length(y))
  }

  # logistic regression on each column of X separately
  for (j in 1:p) {
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
  return(list(betahat = betahat, shat2 = shat2, intercept = intercept))
}

# Compute SER using Wakefields ABF
wakefield_ser <- function(y, betahat, shat2, prior_variance, prior_weights = NULL) {
  # uniform prior on which column of X to select
  p <- length(betahat)
  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / p, p)
  }

  # wakefield ABF
  lbf <- dnorm(betahat, 0, sqrt(prior_variance + shat2), log = TRUE) -
    dnorm(betahat, 0, sqrt(shat2), log = TRUE)
  lbf[is.infinite(shat2)] <- 0 # deal with special case of infinite shat2 (eg happens if X does not vary)

  # normalize
  maxlbf <- max(lbf)
  w <- exp(lbf - maxlbf)

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
    shat2 = shat2
  )
  return(res)
}

# compute SER using GLM
#' @export
fit_glm_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, estimate_prior_variance = F, prior_weights = NULL, family = "binomial") {
  glm_fits <- map_glm(X, y, o, estimate_intercept, family)
  betahat <- glm_fits$betahat
  shat2 <- glm_fits$shat2

  if (estimate_prior_variance) {
    prior_variance <- pmax(0, betahat^2 - shat2)
  }

  res <- wakefield_ser(y, betahat, shat2, prior_variance, prior_weights)
  res$intercept <- glm_fits$intercept
  return(res)
}


shat2_jj <- function(x, y, betahat, intercept) {
  xi <- abs(x * betahat + intercept)
  omega <- pg_kl(rep(1, length(xi)), xi)
  shat2 <- 1 / sum(x^2 * omega)
  return(shat2)
}

#' @export
fit_jjabf_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, estimate_prior_variance = F, prior_weights = NULL) {
  uvb_ser <- fit_uvb_ser(X, y, o,
    prior_variance = 1e10,
    estimate_intercept = estimate_intercept,
    estimate_prior_variance = estimate_prior_variance
  )

  betahat <- uvb_ser$mu
  shat2 <- uvb_ser$var
  intercept <- uvb_ser$intercept

  # replace shat with estimate from JJ approximation
  p <- length(betahat)
  # estimate prior variance
  if (estimate_prior_variance) {
    stop("estimating prior variance not implimented for this SER")
  }

  res <- wakefield_ser(y, betahat, shat2, prior_variance, prior_weights)
  res$intercept <- intercept
  return(res)
}
