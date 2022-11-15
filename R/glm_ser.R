# Adapted from: https://andrewg3311.github.io/susieR_logistic_wflow/susie_logistic_demonstration.html#adjustments_to_susie_code

# compute SER using GLM
fit_glm_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, prior_weights = NULL) {
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
      log.fit <- glm(y ~ X[, j] + 1 + offset(o), family = "binomial") # fit w/ intercept
      intercept[j] <- unname(coef(log.fit)[1])
    } else {
      log.fit <- glm(y ~ X[, j] - 1 + offset(o), family = "binomial") # fit w/out intercept
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


ibss_from_ser <- function(X, y, L = 10, prior_variance = 1., prior_weights = NULL, tol = 1e-3, maxit = 100, estimate_intercept = TRUE, ser_function = NULL) {
  if (is.null(ser_function)) {
    stop("You need to specify a fit function `fit_glm_ser`, `fit_vb_ser`, etc")
  }

  p <- ncol(X)
  n <- nrow(X)

  # place to store posterior info for each l = 1, ..., L
  post_alpha <- matrix(NA, nrow = L, ncol = p)
  post_mu <- matrix(NA, nrow = L, ncol = p)
  post_var <- matrix(NA, nrow = L, ncol = p)
  post_info <- list(alpha = post_alpha, mu = post_mu, var = post_var)

  # store posterior effect estimates
  beta_post_init <- matrix(Inf, nrow = L, ncol = p) # initialize
  beta_post_init2 <- beta_post_init
  beta_post <- matrix(0, nrow = L, ncol = p)

  fixed <- rep(0, n) # fixed portion, estimated from l' != l other SER models

  iter <- 0

  tictoc::tic() # start timer
  fits <- vector(mode = "list", length = L)
  beta_post_history <- vector(mode = "list", length = maxit)
  # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
  while ((norm(beta_post - beta_post_init, "1") > tol) & (norm(beta_post - beta_post_init2, "1") > tol)) {
    beta_post_init2 <- beta_post_init # store from 2 iterations ago
    beta_post_init <- beta_post

    for (l in 1:L) {
      fixed <- fixed - (X %*% beta_post[l, ])[, 1] # remove effect from previous iteration
      ser_l <- ser_function(X, y, o = fixed, prior_variance = prior_variance, estimate_intercept = estimate_intercept, prior_weights = prior_weights)

      # store
      post_info$alpha[l, ] <- ser_l$alpha
      post_info$mu[l, ] <- ser_l$mu
      post_info$var[l, ] <- ser_l$var

      # update beta_post
      beta_post[l, ] <- ser_l$alpha * ser_l$mu

      # add back new fixed portion
      fixed <- fixed + (X %*% beta_post[l, ])[, 1]

      # save current ser
      fits[[l]] <- ser_l
    }
    iter <- iter + 1
    beta_post_history[[iter]] <- beta_post

    if (iter > maxit) {
      warning("Maximum number of iterations reached")
      break
    }
  }
  timer <- tictoc::toc()


  # now, get intercept w/ MLE, holding our final estimate of beta to be fixed
  beta <- colSums(post_info$alpha * post_info$mu)
  pred <- (X %*% beta)[, 1]
  int <- coef(glm(y ~ 1 + offset(pred), family = "binomial"))
  post_info$intercept <- int

  # add the final ser fits
  names(fits) <- paste0("L", 1:L)
  post_info$fits <- fits
  post_info$iter <- iter
  post_info$elapsed_time <- unname(timer$toc - timer$tic)
  post_info$beta_post_history <- beta_post_history
  return(post_info)
}

ibss_monitor_convergence <- function(fit) {
  map_dbl(1:(fit$iter - 1), ~ norm(
    fit$beta_post_history[[.x]] - fit$beta_post_history[[.x + 1]],
    type = "2"
  ))
}
