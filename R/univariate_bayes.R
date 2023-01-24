# Quadrature for "single-effect" logistic regression
# log(p/1-p) = b0 + b1x
# b \sim N(0, s0^2)
# b0 \sim U[-10, 10]

# Misc -------
# sigmoid(x) returns the sigmoid of the elements of x. The sigmoid
# function is also known as the logistic link function. It is the
# inverse of logit(x).
sigmoid <- function(x) {
  y <- x
  y[] <- 0
  y[x > -500] <- 1 / (1 + exp(-x))
  return(y)
}

# logpexp(x) returns log(1 + exp(x)). The computation is performed in a
# Inference -------
# Compute likelihood conditional log likelihood
conditional_log_likelihood <- function(x, y, o, b, b0) {
  psi <- x * b + b0 + o
  ll <- sum(y * psi + logsigmoid(-psi))
  return(ll)
}

# optimize over intercept
null_model_likelihood_with_offset <- function(y, o) {
  x <- runif(length(y))
  f <- function(b0) {
    conditional_log_likelihood(x, y, o, 0, b0)
  }
  return(optimize(f, c(-10, 10), maximum = T)$objective)
}

# 1d quadrature to integrate over effect variable
bayes_logistic_fixed_intercept <- function(x, y, o, b0, s0, quad) {
  # integrate over uniform

  # compute log evidence
  # p(y) = E_p[p(y|b)] = E_q(b)[p(y|b)p(b)/q(b)]  q(b) ~U[-5, 5]
  ll <- purrr::map_dbl(quad$nodes, ~ conditional_log_likelihood(x, y, o, .x, b0))
  logp <- dnorm(quad$nodes, mean = 0, sd = s0, log = T)

  # since we are integrating over uniform need to add back 1/n
  logf <- ll + logp + log(diff(range(quad$nodes)))
  u <- max(logf)
  logpy <- log(sum(exp(logf - u) * quad$weights)) + u

  # compute posterior mean
  # E_b = E_p[p(y|b)] = E_q(b)[p(y|b)p(b)/q(b)]  q(b) ~U[-5, 5]
  mu <- sum(quad$nodes * exp(logf - u) * quad$weights) * exp(u - logpy)
  b2 <- sum(quad$nodes^2 * exp(logf - u) * quad$weights) * exp(u - logpy)

  return(list(logpy = logpy, mu = mu, var = b2 - mu^2))
}

# optionally optimize over intercept
bayes_logistic <- function(x, y, o, b0_init, s0, quad, optimize_b0 = F) {
  if (optimize_b0) {
    # optimize over b0
    f <- function(b0) {
      bayes_logistic_fixed_intercept(x, y, o, b0, s0, quad)$logpy
    }
    b0_opt <- optimise(f, c(-10, 10), maximum = T)$maximum
  } else {
    b0_opt <- b0_init
  }
  # evaluate at optimal fixed intercept b0_opt
  fit <- bayes_logistic_fixed_intercept(x, y, o, b0_opt, s0, quad)
  fit$intercept <- b0_opt
  return(fit)
}

#' Fit logistic SER via quadrature
#'
#' Estimate a fixed intercept using UVB approximation
#' and then numerically integrate over effect variable with normal prior
#' @export
fit_bayes_ser <- function(X, y, o = NULL, prior_variance = 1, estimate_intercept = T, prior_weights = NULL, n = 2^10) {
  # set=up
  p <- dim(X)[2]
  s0 <- sqrt(prior_variance)

  # compute quadrature points once
  quad <- statmod::gauss.quad.prob(n = n, dist = "uniform", l = -10, u = 10)

  # 0 offset if not specified
  if (is.null(o)) {
    o <- rep(0, length(y))
    null_loglik <- sum(dbinom(y, 1, mean(y), log = T))
  } else {
    # fit intercept only model with offsets
    bayes_logistic(X[, 1], y, o = o, 0, 1e-10, quad, optimize_b0 = T)
    null_loglik <- null_model_likelihood_with_offset(y, o)
  }

  # use uvb intercepts as fixed intercept
  uvb <- logisticsusie::fit_uvb_ser(X, y, o = o)
  fits <- bind_rows(purrr::map(1:p, ~ bayes_logistic(X[, .x], y, o, uvb$intercept[.x], s0, quad)))
  fits$lbf <- fits$logpy - null_loglik

  # compute summaries
  alpha <- exp(fits$lbf - matrixStats::logSumExp(fits$lbf))
  lbf_model <- sum(alpha * fits$lbf) - categorical_kl(alpha, rep(1 / p, p))

  # return standard ouput: mu var alpha intercept, lbf
  res <- list(
    mu = fits$mu,
    var = fits$var,
    alpha = alpha,
    intercept = fits$intercept,
    lbf = fits$lbf,
    lbf_model = lbf_model,
    prior_variance = prior_variance,
    logpy = fits$logpy,
    null_loglik = null_loglik
  )
  class(res) <- "ser"
  return(res)
}
