
update_intercept <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + o
  omega <- logisticsusie:::pg_mean(1, xi)
  return(sum(kappa - xb * omega) / sum(omega))
}

update_b <- function(x, y, o, mu, tau, xi, delta, tau0) {
  delta <- delta + o
  omega <- logisticsusie:::pg_mean(1, xi)
  kappa <- y - 0.5
  tau <- sum(omega * x^2) + tau0
  nu <- sum((kappa - omega * delta) * x)
  return(list(mu = nu / tau, tau = tau))
}

update_xi <- function(x, y, o, mu, tau, xi, delta, tau0) {
  delta <- delta + o
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  return(sqrt(xb2))
}

compute_elbo <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  delta <- delta + o
  xb <- (x * mu) + delta
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# include the term that cancels out when xi is up-to-date
compute_elbo2 <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  delta <- delta + o
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi) + 0.5 * omega * (xi^2 - xb2)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# update xi on the fly, so bound is tight
compute_elbo3 <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  delta <- delta + o
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  xi <- sqrt(xb2)
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# deal with offsets so you can compute BFs for each individual SER?
compute_elbo4 <- function(x, y, o, mu, tau, xi, delta, tau0, offset, offset2) {
  kappa <- y - 0.5
  delta <- delta + o
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  xi <- sqrt(xb2)
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

#' Update prior precision of effect
update_tau0 <- function(x, y, o, mu, tau, xi, delta, tau0) {
  b2 <- (mu^2 + 1 / tau)
  return(1 / b2)
}

#' Quick implimentation of univariate VB logistic regression
#' The intercept `delta` is treated as a parameter to be optimized,
#' Normal prior on the effect N(0, 1/tau0)
fit_univariate_vb <- function(x, y, o = 0, delta.init = logodds(mean(y)), tau0 = 1, estimate_intercept = T, maxit = 50, tol = 1e-3) {
  # init
  mu <- 0
  tau <- 1
  delta <- delta.init
  xi <- update_xi(x, y, o, mu, tau, 1, delta, tau0)
  xi <- pmax(xi, 1e-3)

  elbos <- compute_elbo(x, y, o, mu, tau, xi, delta, tau0)
  for (i in seq(maxit)) {
    # rep
    if (estimate_intercept) {
      delta <- update_intercept(x, y, o, mu, tau, xi, delta, tau0)
    }

    b <- update_b(x, y, o, mu, tau, xi, delta, tau0)
    mu <- b$mu
    tau <- b$tau

    xi <- update_xi(x, y, o, mu, tau, xi, delta, tau0)
    elbos <- c(elbos, compute_elbo(x, y, o, mu, tau, xi, delta, tau0))

    if (diff(tail(elbos, 2)) < tol) {
      break
    }
  }

  converged <- diff(tail(elbos, 2)) < tol
  monotone <- logisticsusie:::.monotone(elbos)
  return(list(
    x = x, y = y, o = o,
    mu = mu, tau = tau, xi = xi, delta = delta, tau0 = tau0,
    elbos = elbos,
    converged = converged,
    monotone = monotone
  ))
}

#' Fit a logistic single effect regression model using univariate VB approximation
fit_uvb_ser <- function(X, y, o = NULL, prior_variance = 1.0, intercept.init = logodds(mean(y)), estimate_intercept = T, prior_weights = NULL) {
  tau0 <- 1 / prior_variance
  p <- dim(X)[2]
  if (is.null(o)) {
    # fixed offsets
    o <- 0
  }
  null_model_elbo <- tail(fit_univariate_vb(X[, 1], y, o = o, tau0 = 1e10)$elbos, 1)
  res <- purrr::map(1:p, ~ fit_univariate_vb(X[, .x], y, o = o, tau0 = tau0, delta.init = intercept.init, estimate_intercept = estimate_intercept)) %>%
    dplyr::tibble() %>%
    tidyr::unnest_wider(1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(elbo = tail(elbos, 1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(BF = elbo - null_model_elbo, PIP = exp(elbo - matrixStats::logSumExp(elbo)))
  # return(res)

  lbf_model <- sum(res$BF * res$PIP) - categorical_kl(res$PIP, rep(1 / p, p))
  loglik <- lbf_model + null_model_elbo
  return(list(mu = res$mu, var = 1 / res$tau, alpha = res$PIP, intercept = res$delta, lbf = res$BF, lbf_model = lbf_model, prior_variance = prior_variance, loglik = loglik))
}
