# Updates --------
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

# Single Fit-------------
#' Quick implimentation of univariate VB logistic regression
#' The intercept `delta` is treated as a parameter to be optimized,
#' Normal prior on the effect N(0, 1/tau0)
#' @export
fit_univariate_vb <- function(x, y, o = 0,
                              delta.init = logodds(mean(y) + 1e-10),
                              tau0 = 1,
                              estimate_intercept = T,
                              estimate_prior_variance = F,
                              maxit = 50,
                              tol = 1e-3) {
  # initialize
  mu <- 0
  tau <- 1
  delta <- delta.init
  xi <- update_xi(x, y, o, mu, tau, 1, delta, tau0)
  xi <- pmax(xi, 1e-3)
  elbos <- compute_elbo(x, y, o, mu, tau, xi, delta, tau0)

  # if estimating prior variance, initialize with fixed prior variance fit
  if (estimate_prior_variance) {
    prefit <- fit_univariate_vb(
      x, y, o, delta.init, tau0,
      estimate_intercept = estimate_intercept,
      estimate_prior_variance = F,
      maxit = maxit,
      tol = tol
    )
    mu <- prefit$mu
    tau <- prefit$tau
    delta <- prefit$delta
    xi <- prefit$xi
    elbos <- prefit$elbos
  }

  # loop CAVI update
  for (i in seq(maxit)) {
    if (estimate_intercept) {
      delta <- update_intercept(x, y, o, mu, tau, xi, delta, tau0)
    }

    b <- update_b(x, y, o, mu, tau, xi, delta, tau0)
    mu <- b$mu
    tau <- b$tau

    xi <- update_xi(x, y, o, mu, tau, xi, delta, tau0)
    elbos <- c(elbos, compute_elbo(x, y, o, mu, tau, xi, delta, tau0))

    if (estimate_prior_variance) {
      tau0 <- update_tau0(x, y, o, mu, tau, xi, delta, tau0)
    }

    if (diff(tail(elbos, 2)) < tol) {
      break
    }
  }

  # record convergence
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

# SER--------
#' Fit a logistic single effect regression model using univariate VB approximation
#' @export
fit_uvb_ser <- function(X, y, o = NULL,
                        prior_variance = 1.0,
                        intercept.init = logodds(mean(y) + 1e-10),
                        estimate_intercept = T,
                        estimate_prior_variance = F,
                        prior_weights = NULL) {
  tau0 <- 1 / prior_variance
  p <- dim(X)[2]
  if (is.null(o)) {
    # fixed offsets
    o <- 0
  }

  # null_model_elbo <- tail(fit_univariate_vb(X[, 1], y, o = o, tau0 = 1e10)$elbos, 1)
  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))

  res <- purrr::map(1:p, ~ fit_univariate_vb(
    X[, .x], y,
    o = o,
    tau0 = tau0,
    delta.init = intercept.init,
    estimate_intercept = estimate_intercept,
    estimate_prior_variance = estimate_prior_variance
  )) %>%
    dplyr::tibble() %>%
    tidyr::unnest_wider(1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(elbo = tail(elbos, 1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      lbf = elbo - null_likelihood,
      alpha = exp(lbf - matrixStats::logSumExp(lbf))
    )

  lbf_model <- sum(res$lbf * res$alpha) - categorical_kl(res$alpha, rep(1 / p, p))
  loglik <- lbf_model + null_likelihood
  res <- list(
    mu = res$mu,
    var = 1 / res$tau,
    alpha = res$alpha,
    intercept = res$delta,
    lbf = res$lbf,
    elbo = res$elbo,
    lbf_model = lbf_model,
    prior_variance = 1 / res$tau0,
    loglik = loglik,
    null_loglik = null_likelihood,
    o = o
  )
  return(res)
}
