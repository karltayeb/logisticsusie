# Just a simple wrapper to binsusie for fitting the logistic SER
# we do this to have a common interface with the `vb_ser` and `glm_ser`

#' Compute log Bayes Factors for each feature
#' Takes a SER and give back the log-Bayes Factors for each feature
compute_lbf_variable <- function(ser, l = 1) {
  p <- dim(ser$data$X)[2]
  elbo <- with(ser, purrr::map_dbl(1:p, ~ compute_elbo2(
    x = data$X[, .x],
    y = data$y, o = 0,
    mu = params$mu[l, .x],
    tau = 1 / params$var[l, .x],
    xi = params$xi,
    delta = params$delta[1, 1],
    tau0 = 1 / hypers$prior_variance
  )))
  null_model_elbo <- tail(fit_univariate_vb(ser$data$X[, 1], ser$data$y, tau0 = 1e10)$elbos, 1)
  return(elbo - null_model_elbo)
}

#' Fit a logistic single effect regression with our variational approximation
#' Returns a list with posterior summary
fit_bin_ser <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, prior_weights = NULL) {
  fit <- binsusie(X, y, L = 1, center = F, prior_variance = prior_variance, estimate_prior_variance = F)
  lbf <- compute_lbf_variable(fit)
  # also just elbo - null_model_elbo
  lbf_model <- sum(lbf * fit$pip) - categorical_kl(fit$pip, rep(1 / 20, 20))
  loglik <- tail(fit$elbo, 1)
  intercept <- rep(fit$params$delta[1, 1], length(lbf))
  return(list(
    alpha = fit$alpha, mu = fit$mu, var = fit$params$var, intercept = intercept,
    lbf = lbf, lbf_model = lbf_model,
    prior_variance = prior_variance, loglik = loglik
  ))
}


#' reestimate xi for each feature, this is the "corrected" SER
compute_corrected_lbf <- function(X, y, fit) {
  p <- dim(X)[2]
  elbo <- purrr::map_dbl(1:p, ~ compute_elbo3(
    x = X[, .x],
    y = y, o = 0,
    mu = fit$mu[.x],
    tau = 1 / fit$var[.x],
    xi = 0,
    delta = fit$intercept[.x],
    tau0 = 1 / fit$prior_variance
  ))
  null_model_elbo <- tail(fit_univariate_vb(X[, 1], y, tau0 = 1e10)$elbos, 1)
  return(elbo - null_model_elbo)
}

#' Correct a binomial SER fit by recomputing the lbfs
correct_bin_ser <- function(X, y, fit) {
  # compute log BFs using optimal xi for each feature
  lbf <- compute_corrected_lbf(X, y, fit)

  # update alpha with new approximate of Bayes Factors
  alpha <- exp(lbf - matrixStats::logSumExp(lbf))

  # compute logBF for the entire SER model
  p <- length(lbf)
  lbf_model <- sum(lbf * fit$alpha) - categorical_kl(fit$alpha, rep(1 / p, p))

  # use logBF + null_likelihood to compute log likelihood
  null_model_elbo <- tail(fit_univariate_vb(X[, 1], y, tau0 = 1e10)$elbos, 1)
  loglik <- lbf_model + null_model_elbo

  # update
  fit$alpha <- alpha
  fit$lbf <- lbf
  fit$lbf_model <- lbf_model
  fit$loglik <- loglik
  return(fit)
}

#' Fit a logistic single effect regression with our variational approximation
#' Returns a list with posterior summary
fit_bin_ser_corrected <- function(X, y, o = NULL, prior_variance = 1.0, estimate_intercept = T, prior_weights = NULL) {
  fit <- fit_bin_ser(X, y, o, prior_variance, estimate_intercept, prior_weigts)
  fit <- correct_bin_ser(X, y, fit)
}
