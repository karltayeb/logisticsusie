
update_intercept <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu)
  omega <- logisticsusie:::pg_mean(1, xi)
  return(sum(kappa - xb * omega) / sum(omega))
}

update_b <- function(x, y, mu, tau, xi, delta, tau0) {
  omega <- logisticsusie:::pg_mean(1, xi)
  kappa <- y - 0.5
  tau <- sum(omega * x^2) + tau0
  nu <- sum((kappa - omega * delta) * x)
  return(list(mu = nu / tau, tau = tau))
}

update_xi <- function(x, y, mu, tau, xi, delta, tau0) {
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  return(sqrt(xb2))
}

compute_elbo <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + delta
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# include the term that cancels out when xi is up-to-date
compute_elbo2 <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi) + 0.5 * omega * (xi^2 - xb2)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# update xi on the fly, so bound is tight
compute_elbo3 <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  xi <- sqrt(xb2)
  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * xb) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# compute the elbo explicitly
compute_elbo3 <- function(x, y, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + delta
  xb2 <- x^2 * (mu^2 + 1 / tau) + 2 * x * mu * delta + delta^2
  omega <- logisticsusie:::pg_mean(1, xi)

  elbo <- sum(-log(2) + kappa * xb - 0.5 * xb2 * omega) -
    logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0) -
    sum(logisticsusie:::pg_kl(1, xi))
  return(elbo)
}

#' Update prior precision of effect
update_tau0 <- function(x, y, mu, tau, xi, delta, tau0) {
  b2 <- (mu^2 + 1 / tau)
  return(1 / b2)
}

#' Quick implimentation of univariate VB logistic regression
#' The intercept `delta` is treated as a parameter to be optimized,
#' Normal prior on the effect N(0, 1/tau0)
fit_univariate_vb <- function(x, y, delta.init = 0, tau0 = 1, estimate_intercept = T, maxit = 50, tol = 1e-3) {
  # init
  mu <- 0
  tau <- 1
  delta <- delta.init
  xi <- update_xi(x, y, mu, tau, 1, delta, tau0)
  xi <- pmax(xi, 1e-3)

  elbos <- compute_elbo(x, y, mu, tau, xi, delta, tau0)
  for (i in seq(maxit)) {
    # rep
    if (estimate_intercept) {
      delta <- update_intercept(x, y, mu, tau, xi, delta, tau0)
    }

    b <- update_b(x, y, mu, tau, xi, delta, tau0)
    mu <- b$mu
    tau <- b$tau

    xi <- update_xi(x, y, mu, tau, xi, delta, tau0)
    elbos <- c(elbos, compute_elbo(x, y, mu, tau, xi, delta, tau0))

    if (diff(tail(elbos, 2)) < tol) {
      break
    }
  }

  converged <- diff(tail(elbos, 2)) < tol
  monotone <- logisticsusie:::.monotone(elbos)
  return(list(
    x = x, y = y,
    mu = mu, tau = tau, xi = xi, delta = delta, tau0 = tau0,
    elbos = elbos,
    converged = converged,
    monotone = monotone
  ))
}
