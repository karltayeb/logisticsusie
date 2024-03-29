# Univariate VB Updates --------
update_intercept_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + o
  omega <- pg_mean(1, xi)
  return(sum(kappa - xb * omega) / sum(omega))
}

update_b_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  delta <- delta + o
  omega <- pg_mean(1, xi)
  kappa <- y - 0.5
  tau <- sum(omega * x^2) + tau0
  nu <- sum((kappa - omega * delta) * x)
  return(list(mu = nu / tau, tau = tau))
}

compute_psi_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  psi <- x * mu + delta + o
  return(psi)
}

compute_psi2_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  Vo <- o2 - o^2
  Vxb <- x^2 / tau

  xb <- x * mu
  psi2 <- (xb + delta + o)^2 + Vo + Vxb

  if(any(is.nan(sqrt(psi2)))){browser()}
  return(psi2)
}

update_xi_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  psi2 <- compute_psi2_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
  return(sqrt(psi2))
}

compute_elbo_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  psi <- compute_psi_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi)
  kl <- normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# include the term that cancels out when xi is up-to-date
compute_elbo_vb2 <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  delta <- delta + o
  psi <- compute_psi_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
  psi2 <- compute_psi2_vb(x, y, o, o2, mu, tau, xi, delta, tau0)

  omega <- pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi) + 0.5 * omega * (xi^2 - psi2)
  kl <- normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# update xi on the fly, so bound is tight
compute_elbo_vb3 <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  psi <- compute_psi_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
  psi2 <- compute_psi2_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
  xi <- sqrt(psi2)
  omega <- pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi)
  kl <- normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

#' Update prior precision of effect
update_tau0_vb <- function(x, y, o, o2, mu, tau, xi, delta, tau0) {
  b2 <- (mu^2 + 1 / tau)
  return(1 / b2)
}

# Single Fit-------------
#' Quick implimentation of univariate VB logistic regression
#' The intercept `delta` is treated as a parameter to be optimized,
#' Normal prior on the effect N(0, 1/tau0)
#' @export
fit_univariate_vb <- function(x, y, o = 0, o2 = 0,
                              mu.init = 0,
                              tau.init = 1,
                              xi.init = NULL,
                              delta.init = logodds(mean(y) + 1e-10),
                              tau0 = 1,
                              estimate_intercept = T,
                              maxit = 50,
                              tol = 1e-3) {
  # init
  mu <- mu.init
  tau <- tau.init
  delta <- delta.init
  if(is.null(xi.init)){
    xi <- update_xi_vb(x, y, o, o2, mu, tau, 1, delta, tau0)
    xi <- pmax(xi, 1e-3)
  }

  elbos <- compute_elbo_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
  for (i in seq(maxit)) {
    # rep
    if (estimate_intercept) {
      delta <- update_intercept_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
    }

    b <- update_b_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
    mu <- b$mu
    tau <- b$tau

    xi <- update_xi_vb(x, y, o, o2, mu, tau, xi, delta, tau0)
    elbos <- c(elbos, compute_elbo_vb(x, y, o, o2, mu, tau, xi, delta, tau0))

    if (diff(tail(elbos, 2)) < tol) {
      break
    }
  }

  converged <- diff(tail(elbos, 2)) < tol
  monotone <- .monotone(elbos)
  return(list(
    x = x, y = y, o = o, o2=o2,
    mu = mu, tau = tau, xi = xi, delta = delta, tau0 = tau0,
    elbos = elbos,
    converged = converged,
    monotone = monotone
  ))
}

# SER--------
#' Fit a logistic single effect regression model using univariate VB approximation
#' @export
fit_uvb_ser <- function(X, y, o = NULL, o2=NULL,
                        prior_variance = 1.0,
                        intercept.init = logodds(mean(y) + 1e-10),
                        estimate_intercept = T,
                        prior_weights = NULL,
                        mapper = purrr::map) {
  tau0 <- 1 / prior_variance
  p <- dim(X)[2]
  if (is.null(o)) {
    # fixed offsets
    o <- 0
  }
  if (is.null(o2)) {
    o2 <- o^2
  }

  #null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  null_likelihood <- tail(fit_univariate_vb(X[, 1], y, o = o, o2=o2, tau0 = 1e10)$elbos, 1)

  res <- mapper(1:p, ~ fit_univariate_vb(
    X[, .x], y,
    o = o, o2 = o2,
    tau0 = tau0,
    delta.init = intercept.init,
    estimate_intercept = estimate_intercept
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
    lbf_model = lbf_model,
    prior_variance = prior_variance,
    loglik = loglik,
    null_loglik = null_likelihood
  )
  return(res)
}
