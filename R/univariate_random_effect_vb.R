# Modify the univariate logistic VB to include a random effect offset
# instead of providing a fixed offset, provide first and second moments of the offset term

# Updates-------
update_intercept_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  xb <- (x * mu) + o$mu
  omega <- logisticsusie:::pg_mean(1, xi)
  return(sum(kappa - xb * omega) / sum(omega))
}

update_b_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  delta <- delta + o$mu
  omega <- logisticsusie:::pg_mean(1, xi)
  kappa <- y - 0.5
  tau <- sum(omega * x^2) + tau0 # posterior precision
  nu <- sum((kappa - omega * delta) * x) # unormalized posterior mean
  return(list(mu = nu / tau, tau = tau))
}

compute_psi2_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  Vo <- o$mu2 - o$mu^2
  Vxb <- x^2 * (1 / tau)
  xb <- x * mu
  o <- o$mu
  psi2 <- (xb + delta + o)^2 + Vo + Vxb
  return(psi2)
}

update_xi_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  psi2 <- compute_psi2_re(x, y, o, mu, tau, xi, delta, tau0)
  return(sqrt(psi2))
}

compute_elbo_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  delta <- delta + o$mu
  psi <- (x * mu) + delta
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# include the term that cancels out when xi is up-to-date
compute_elbo2_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  delta2 <- o$mu2 + delta^2 + (2 * o$mu * delta)
  delta <- delta + o$mu

  xb <- (x * mu)
  psi <- xb + delta

  xb2 <- x^2 * (mu^2 + 1 / tau)
  psi2 <- xb2 + delta2 + 2 * xb * delta

  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi) + 0.5 * omega * (xi^2 - psi2)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}


# Single fit --------
#' Quick implimentation of univariate VB logistic regression
#' The intercept `delta` is treated as a parameter to be optimized,
#' Normal prior on the effect N(0, 1/tau0)
fit_univariate_vb_re <- function(x, y, o = NULL, delta.init = logodds(mean(y) + 1e-10), tau0 = 1, estimate_intercept = T, estimate_prior_variance = F, maxit = 50, tol = 1e-3) {
  # if no random offset is provided, same as univariate vb
  if (is.null(o)) {
    o <- list(mu = 0, mu2 = 0)
  }

  mu <- 0
  tau <- 1
  delta <- delta.init
  xi <- update_xi_re(x, y, o, mu, tau, 1, delta, tau0)
  xi <- pmax(xi, 1e-3)
  elbos <- compute_elbo_re(x, y, o, mu, tau, xi, delta, tau0)

  # if estimating prior variance, initialize with fixed prior variance fit
  if (estimate_prior_variance) {
    prefit <- fit_univariate_vb_re(
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

  # loop CAVI
  for (i in seq(maxit)) {
    # rep
    if (estimate_intercept) {
      delta <- update_intercept_re(x, y, o, mu, tau, xi, delta, tau0)
    }

    b <- update_b_re(x, y, o, mu, tau, xi, delta, tau0)
    mu <- b$mu
    tau <- b$tau

    xi <- update_xi_re(x, y, o, mu, tau, xi, delta, tau0)
    elbos <- c(elbos, compute_elbo_re(x, y, o, mu, tau, xi, delta, tau0))

    if (estimate_prior_variance) {
      tau0 <- update_tau0(x, y, o, mu, tau, xi, delta, tau0)
    }

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

# SER--------

compute_psi2_ser_re <- function(X, alpha, mu, var, delta) {
  # expected predictions
  Xb <- Matrix::drop(X %*% (alpha * mu))
  Zd <- sum(alpha * delta)
  psi <- Xb + Zd

  # variance components
  b2 <- alpha * (mu^2 + var)
  Xb2 <- Matrix::drop(X^2 %*% b2)
  VXb <- Xb2 - Xb^2

  # second moment
  psi2 <- psi^2 + VXb
  return(psi2)
}

#' Fit a logistic single effect regression model using univariate VB approximation
#' @export
fit_uvb_ser_re <- function(X, y, o = NULL, prior_variance = 1.0, intercept.init = logodds(mean(y) + 1e-10), estimate_intercept = T, estimate_prior_variance = F, prior_weights = NULL) {
  tau0 <- 1 / prior_variance
  p <- dim(X)[2]
  if (is.null(o)) {
    # fixed offsets
    o <- list(mu = 0, mu2 = 0)
  }
  null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  res <- purrr::map(1:p, ~ fit_univariate_vb_re(X[, .x], y, o = o, tau0 = tau0, delta.init = intercept.init, estimate_intercept = estimate_intercept, estimate_prior_variance = estimate_prior_variance)) %>%
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

  psi <- drop(X %*% (res$alpha * res$mu)) + sum(res$alpha * res$delta)
  psi2 <- compute_psi2_ser_re(X, res$alpha, res$mu, 1 / res$tau, res$delta)

  res <- list(
    mu = res$mu,
    var = 1 / res$tau,
    alpha = res$alpha,
    intercept = res$delta,
    lbf = res$lbf,
    elbo = res$elbo,
    lbf_model = lbf_model,
    prior_variance = prior_variance,
    loglik = loglik,
    null_loglik = null_likelihood,
    o = o,
    psi = psi,
    psi2 = psi2
  )
  return(res)
}



# IBSS2 ---------

#' Add the predictions from an SER to the random offset
add_re <- function(o, ser, X) {
  # add to offset
  o$mu2 <- ser$psi2 + o$mu2 + 2 * o$mu * ser$psi
  o$mu <- o$mu + ser$psi
  return(o)
}

#' Subtract the predictions from an SER to the random offset
sub_re <- function(o, ser, X) {
  # remove from offset
  o$mu <- o$mu - ser$psi
  o$mu2 <- o$mu2 - ser$psi2 - 2 * o$mu * ser$psi
  return(o)
}

valid_offset <- function(o) {
  return(all(o$mu2 - o$mu^2 >= 0))
}


#' 2 moments Iterative Bayesian Stepwise Selectioin
#' @export
ibss2m <- function(X, y, L = 10, prior_variance = 1., prior_weights = NULL, tol = 1e-3, maxit = 100, estimate_intercept = TRUE) {
  # data dimensions
  p <- ncol(X)
  n <- nrow(X)

  # place to store posterior info for each l = 1, ..., L
  alpha <- matrix(NA, nrow = L, ncol = p)
  mu <- matrix(NA, nrow = L, ncol = p)
  var <- matrix(NA, nrow = L, ncol = p)

  # store posterior effect estimates
  beta_post <- matrix(0, nrow = L, ncol = p)

  # fixed portion, estimated from l' != l other SER models
  fixed <- list(mu = rep(0, n), mu2 = rep(0, n))
  iter <- 0

  tictoc::tic() # start timer

  # initialize empty SERs
  empty_ser <- list(
    alpha = rep(1 / p, p),
    mu = rep(0, p),
    var = (rep(0, p)),
    intercept = 0,
    psi = rep(0, n),
    psi2 = rep(0, n)
  )

  fits <- vector(mode = "list", length = L)
  for (l in 1:L) {
    fits[[l]] <- empty_ser
  }
  beta_post_history <- vector(mode = "list", length = maxit)
  # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
  while (ibss_l2(beta_post_history, iter - 1) > tol) {
    for (l in 1:L) {
      # remove effect from previous iteration
      ser_l <- fits[[l]]
      fixed <- sub_re(fixed, ser_l, X)

      # fit new SER with random offset
      ser_l <- fit_uvb_ser_re(X, y,
        o = fixed,
        prior_variance = prior_variance,
        estimate_intercept = estimate_intercept,
        prior_weights = prior_weights
      )

      # store
      alpha[l, ] <- ser_l$alpha
      mu[l, ] <- ser_l$mu
      var[l, ] <- ser_l$var

      # update beta_post
      beta_post[l, ] <- ser_l$alpha * ser_l$mu

      # add back new (random) fixed portion
      fixed <- add_re(fixed, ser_l, X)

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
  beta <- colSums(alpha * mu)
  pred <- (X %*% beta)[, 1]
  int <- coef(glm(y ~ 1 + offset(pred), family = "binomial"))

  # add the final ser fits
  names(fits) <- paste0("L", 1:L)

  res <- list(
    alpha = alpha,
    mu = mu,
    var = var,
    intercept = int,
    fits = fits,
    iter = iter,
    elapsed_time = unname(timer$toc - timer$tic),
    beta_post_history = head(beta_post_history, iter),
    psi = fixed
  )
  return(res)
}
