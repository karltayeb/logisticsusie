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

compute_psi_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  xb <- x * mu
  psi <- xb + delta + o$mu
  return(psi)
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
  psi <- compute_psi_re(x, y, o, mu, tau, xi, delta, tau0)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# include the term that cancels out when xi is up-to-date
compute_elbo2_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  psi <- compute_psi_re(x, y, o, mu, tau, xi, delta, tau0)
  psi2 <- compute_psi2_re(x, y, o, mu, tau, xi, delta, tau0)

  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi) + 0.5 * omega * (xi^2 - psi2)
  kl <- logisticsusie:::normal_kl(mu, 1 / tau, 0, 1 / tau0)
  return(sum(bound) - kl)
}

# update xi on-the-fly
compute_elbo3_re <- function(x, y, o, mu, tau, xi, delta, tau0) {
  kappa <- y - 0.5
  psi <- compute_psi_re(x, y, o, mu, tau, xi, delta, tau0)
  psi2 <- compute_psi2_re(x, y, o, mu, tau, xi, delta, tau0)
  xi <- sqrt(psi2) # optimal xi

  omega <- logisticsusie:::pg_mean(1, xi)
  bound <- log(sigmoid(xi)) + (kappa * psi) - (0.5 * xi) #+ 0.5 * omega * (xi^2 - psi2)
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
    if (estimate_prior_variance) {
      tau0 <- update_tau0(x, y, o, mu, tau, xi, delta, tau0)
    }

    elbos <- c(elbos, compute_elbo_re(x, y, o, mu, tau, xi, delta, tau0))
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
  # NOTE: DOES NOT INCLUDE OFFSET!
  psi2 <- psi^2 + VXb
  return(psi2)
}

#' Fit a logistic single effect regression model using univariate VB approximation
#' @export
fit_uvb_ser_re <- function(X, y, o = NULL,
                           prior_variance = 1.0,
                           intercept.init = logodds(mean(y) + 1e-10),
                           estimate_intercept = T,
                           estimate_prior_variance = F,
                           prior_weights = NULL) {
  tau0 <- 1 / prior_variance
  X <- as.matrix(X)
  p <- dim(X)[2]
  if (is.null(o)) {
    # fixed offsets
    o <- list(mu = 0, mu2 = 0)
    null_likelihood <- sum(dbinom(y, 1, mean(y), log = T))
  } else {
    # if there is a random effect fit the null model accounting for that
    null_model <- fit_univariate_vb_re(X[, 1], y,
      o = o, tau0 = 1e10,
      delta.init = intercept.init,
      estimate_intercept = estimate_intercept,
      estimate_prior_variance = F
    )
    null_likelihood <- tail(null_model$elbos, 1)
  }

  # fit uvb_re for each column of X
  res <- purrr::map(1:p, ~ fit_univariate_vb_re(
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

  # compute LBF for model
  lbf_model <- sum(res$lbf * res$alpha) -
    categorical_kl(res$alpha, rep(1 / p, p))
  loglik <- lbf_model + null_likelihood

  psi <- drop(X %*% (res$alpha * res$mu)) + sum(res$alpha * res$delta)
  psi2 <- compute_psi2_ser_re(X, res$alpha, res$mu, 1 / res$tau, res$delta)

  res2 <- list(
    mu = res$mu,
    var = 1 / res$tau,
    alpha = res$alpha,
    intercept = res$delta,
    lbf = res$lbf,
    elbo = res$elbo,
    lbf_model = lbf_model,
    elbo_model = lbf_model + null_likelihood,
    prior_variance = ifelse(estimate_prior_variance,
      list(1 / res$tau0),
      list(prior_variance)
    )[[1]],
    loglik = loglik,
    null_loglik = null_likelihood,
    o = o,
    psi = psi,
    psi2 = psi2
  )
  return(res2)
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


jj_bound_re <- function(X, y, psi) {
  xi <- sqrt(psi$mu2)
  Xb <- psi$mu
  kappa <- y - 0.5
  n <- 1
  bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * xi)
  return(bound)
}

compute_elbo_susie_re <- function(X, y, psi, fits, prior_weights) {
  jj <- jj_bound_re(X, y, psi) # E[log p(y, w | b) - logp(w)]
  kl <- purrr::map_dbl(fits, ~ .compute_ser_kl(
    .x$alpha, prior_weights, .x$mu, .x$var, 0, .x$prior_var
  )) # KL(q(b) || p(b))
  elbo <- sum(jj) - sum(kl)
  return(elbo)
}

# GLOBAL AND SER-centric BOUNDS

# Compute JJ bound using global optimum of xi
ibss2_compute_jj <- function(y, post) {
  # compute overall prediction moments
  psi <- colSums(post$psi)
  psi2 <- colSums(post$psi)^2 + colSums(post$psi2 - post$psi^2)

  # optimal xi
  xi <- sqrt(psi2)

  # global jj bound
  jj <- log(sigmoid(xi)) - 0.5 * xi + (y - 0.5) * psi
  return(sum(jj))
}

# compute optimal jj bound partitioning over \gamma_l
ibss2_compute_jj_l <- function(X, y, post, l) {
  n <- dim(X)[1]
  psi_ml <- colSums(post$psi[-l, ])
  Vpsi_ml <- colSums(post$psi2[-l, ] - post$psi[-l, ]^2)

  psi <- do.call(rbind, purrr::map(
    1:n, ~ X[.x, ] * post$mu[l, ] + post$delta[l, ] + psi_ml[.x]
  ))
  Vpsi <- do.call(rbind, purrr::map(
    1:n, ~ X[.x, ]^2 * post$var[l, ] + Vpsi_ml[.x]
  ))
  psi2 <- psi^2 + Vpsi

  xi <- sqrt(psi2)
  jj <- log(sigmoid(xi)) - 0.5 * xi + (y - 0.5) * psi
  jj <- drop(jj %*% post$alpha[l, ])
  return(sum(jj))
}

#' compute KL divergence of SERs
ibss2_compute_kl <- function(post) {
  L <- dim(post$mu)[1]
  kl <- purrr::map_dbl(1:L, ~ .compute_ser_kl(
    alpha = post$alpha[.x, ],
    mu = post$mu[.x, ],
    var = post$var[.x, ],
    pi = post$pi[.x, ],
    prior_mean = 0, var0 = post$prior_variance[.x, ]
  ))
  return(kl)
}

#' compute ELBO for ibss2 fit
ibss2_compute_elbo <- function(y, post) {
  sum(ibss2_compute_jj(y, post)) - sum(ibss2_compute_kl(post))
}

#' compute global and SER-centric ELBOs
compute_all_elbos <- function(X, y, post) {
  L <- dim(post$mu)[1]
  jj_global <- ibss2_compute_jj(y, post)
  jj_l <- purrr::map_dbl(1:L, ~ ibss2_compute_jj_l(X, y, post, .x))
  kl <- sum(ibss2_compute_kl(post))
  elbos <- c(jj_global, jj_l) - kl
  names(elbos) <- c("global", paste0("L", 1:L))
  return(list(elbos))
}

#' extract posterior summaries from list of uvb_re fits
make_post_re <- function(fits) {
  p <- length(fits[[1]]$mu)
  L <- length(fits)
  post <- list(
    mu = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "mu"))),
    var = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "var"))),
    alpha = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "alpha"))),
    pi = matrix(rep(1 / p, L * p), nrow = L),
    prior_variance = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "prior_variance"))),
    delta = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "intercept"))),
    psi = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "psi"))),
    psi2 = do.call(rbind, purrr::map(fits, ~ purrr::pluck(.x, "psi2")))
  )
  return(post)
}

#' 2 moments Iterative Bayesian Stepwise Selectioin
#' @export
ibss2m <- function(X, y, L = 10, prior_variance = 1., prior_weights = NULL, tol = 1e-3, maxit = 100, estimate_intercept = T, estimate_prior_variance = F) {
  # data dimensions
  p <- ncol(X)
  n <- nrow(X)

  if (is.null(prior_weights)) {
    prior_weights <- rep(1 / p, p)
  }

  # place to store posterior info for each l = 1, ..., L
  alpha <- matrix(NA, nrow = L, ncol = p)
  mu <- matrix(NA, nrow = L, ncol = p)
  var <- matrix(NA, nrow = L, ncol = p)

  # store posterior effect estimates
  elbos <- c(-Inf)

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

  # initialize
  fits <- vector(mode = "list", length = L)
  for (l in 1:L) {
    fits[[l]] <- empty_ser
  }

  # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
  for (iter in 1:maxit) {
    for (l in 1:L) {
      # remove effect from previous iteration
      ser_l <- fits[[l]]
      fixed <- sub_re(fixed, ser_l, X)

      # fit new SER with random offset
      ser_l <- fit_uvb_ser_re(X, y,
        o = fixed,
        prior_variance = prior_variance,
        estimate_intercept = estimate_intercept,
        estimate_prior_variance = estimate_prior_variance,
        prior_weights = prior_weights
      )

      # add back new (random) fixed portion
      fixed <- add_re(fixed, ser_l, X)

      # save current ser
      fits[[l]] <- ser_l
    }

    post <- make_post_re(fits)
    elbos <- c(elbos, compute_all_elbos(X, y, post))

    # call converged if all bounds converged
    # NOTE: updates are not necessarily monotonic
    # this is just a measure of how much q has changed since last iteration
    elbo_diff <- tail(elbos, 2)
    elbo_diff <- max(elbo_diff[[2]] - elbo_diff[[1]])
    if (elbo_diff < tol) {
      break
    }

    if (iter > maxit) {
      warning("Maximum number of iterations reached")
      break
    }
  }
  timer <- tictoc::toc()


  # now, get intercept w/ MLE, holding our final estimate of beta to be fixed
  beta <- colSums(post$alpha * post$mu)
  pred <- (X %*% beta)[, 1]
  int <- coef(glm(y ~ 1 + offset(pred), family = "binomial"))

  # add the final ser fits
  names(fits) <- paste0("L", 1:L)

  # post$pip <- get_pip(post$alpha)
  # post$elbos <- as.data.frame(do.call(rbind, tail(elbos, -1)))
  # post$elapsed_time <- unname(timer$toc - timer$tic)
  # return(post)

  elbos <- as.data.frame(do.call(rbind, tail(elbos, -1)))
  res <- list(
    alpha = post$alpha,
    mu = post$mu,
    var = post$var,
    pip = get_pip(post$alpha),
    intercept = int,
    fits = fits,
    iter = iter,
    elapsed_time = unname(timer$toc - timer$tic),
    elbo = elbos,
    psi = fixed,
    converged = iter < maxit,
    # monotone = .monotone(elbos),
    tol = tol
  )
  return(res)
}
