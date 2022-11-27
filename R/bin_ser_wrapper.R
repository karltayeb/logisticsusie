# Just a simple wrapper to binsusie for fitting the logistic SER
# we do this to have a common interface with the `vb_ser` and `glm_ser`


compute_elbo_idx <- function(ser, o = 0, idx) {
  return(with(ser, compute_elbo2(
    x = data$X[, idx],
    y = data$y, o = o,
    mu = params$mu[idx],
    tau = 1 / params$var[idx],
    xi = params$xi,
    delta = params$delta,
    tau0 = 1 / hypers$prior_variance
  )))
}

# Binomial SER ---------
#' Compute log Bayes Factors for each feature
#' Takes a SER and give back the log-Bayes Factors for each feature
#' @export
compute_lbf_variable <- function(X, y, o, fit) {
  p <- dim(X)[2]
  intercept <- fit$intercept
  if (length(intercept) == 1) {
    intercept <- rep(intercept, p)
  }
  elbo <- purrr::map_dbl(1:p, ~ compute_elbo2(
    x = X[, .x],
    y = y, o = o,
    mu = fit$mu[.x],
    tau = 1 / fit$var[.x],
    xi = fit$xi,
    delta = intercept[.x],
    tau0 = 1 / fit$prior_variance
  ))
  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  return(elbo - null_likelihood)
}

#' Fit a logistic single effect regression with our variational approximation
#' Returns a list with posterior summary
#' @export
fit_bin_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, estimate_prior_variance = F, prior_weights = NULL) {
  if (is.null(o)) {
    o <- 0
  }
  fit <- binser(X, y, o = o, center = F, prior_variance = prior_variance, estimate_prior_variance = estimate_prior_variance)
  p <- dim(X)[2]
  loglik <- tail(fit$elbo, 1)
  intercept <- rep(fit$params$delta, )
  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  fit <- list(
    alpha = fit$params$alpha,
    mu = fit$params$mu,
    var = fit$params$var,
    intercept = intercept,
    xi = fit$params$xi,
    prior_variance = fit$hypers$prior_variance,
    loglik = loglik, null_loglik = null_likelihood
  )

  lbf <- compute_lbf_variable(X, y, o, fit)
  lbf_model <- loglik - null_likelihood

  fit$lbf <- lbf
  fit$lbf_model <- lbf_model
  return(fit)
}


# Binomial SER + correction -----------
#' reestimate xi for each feature, this is the "corrected" SER
#' @export
compute_corrected_lbf <- function(X, y, o, fit) {
  p <- dim(X)[2]
  intercept <- fit$intercept
  if (length(intercept) == 1) {
    inercept <- rep(intercept, p)
  }
  elbo <- purrr::map_dbl(1:p, ~ compute_elbo3(
    x = X[, .x],
    y = y, o = o,
    mu = fit$mu[.x],
    tau = 1 / fit$var[.x],
    xi = 0,
    delta = inercept[.x],
    tau0 = 1 / fit$prior_variance
  ))

  # null_model_elbo <- tail(fit_univariate_vb(X[, 1], y, tau0 = 1e10)$elbos, 1)
  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  return(elbo - null_likelihood)
}

#' Correct a binomial SER fit by recomputing the lbfs
#' @export
correct_bin_ser <- function(X, y, o, fit) {
  # compute log BFs using optimal xi for each feature
  lbf <- compute_corrected_lbf(X, y, o, fit)

  # update alpha with new approximate of Bayes Factors
  alpha <- exp(lbf - matrixStats::logSumExp(lbf))

  # compute logBF for the entire SER model
  p <- length(lbf)
  lbf_model <- sum(lbf * alpha) - categorical_kl(alpha, rep(1 / p, p))

  # use logBF + null_likelihood to compute log likelihood
  # null_model_elbo <- tail(fit_univariate_vb(X[, 1], y, tau0 = 1e10)$elbos, 1)
  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  loglik <- lbf_model + null_likelihood

  # update
  fit$alpha <- alpha
  fit$lbf <- lbf
  fit$lbf_model <- lbf_model
  fit$loglik <- loglik
  return(fit)
}

#' Fit a logistic single effect regression with our variational approximation
#' Returns a list with posterior summary
#' @export
fit_bin_ser_corrected <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, prior_weights = NULL) {
  if (is.null(o)) {
    o <- 0
  }
  fit <- fit_bin_ser(X, y, o, prior_variance, estimate_intercept, prior_weigts)
  fit <- correct_bin_ser(X, y, o, fit)
  return(fit)
}



# VEB Boost SER --------
# Note: if we need more from this figure out
#   - how to incorporate offset
#   - and fix the prior variance
#' @export
fit_veb_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, prior_weights = NULL) {
  if (!is.null(o)) {
    stop("offset not implimented for VEB-boost SER")
  }

  veb.fit <- VEB.Boost::veb_boost_stumps(
    X, y,
    family = "binomial",
    include_stumps = FALSE,
    max_log_prior_var = 35,
    scale_X = "NA",
    growMode = "NA",
    changeToConstant = F,
    k = 1
  )

  alpha <- t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$learner$currentFit$alpha)))[1, ]
  mu <- t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$learner$currentFit$mu)))[1, ]
  var <- t(do.call(cbind, lapply(veb.fit$leaves, function(x) x$learner$currentFit$sigma2_post)))[1, ]
  elbo <- tail(veb.fit$ELBO_progress[[2]], 1)

  lbf <- log(alpha)
  p <- length(alpha)
  c <- elbo - (sum(alpha * log(alpha)) - categorical_kl(alpha, rep(1 / p, p)))
  lbf <- log(alpha) + c

  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  lbf_model <- elbo - null_likelihood

  res <- list(
    alpha = alpha,
    mu = mu,
    var = var,
    lbf = lbf,
    lbf_model = lbf_model,
    prior_variance = veb.fit$learner$currentFit$V,
    loglik = elbo,
    null_loglik = null_likelihood
  )
  return(res)
}
