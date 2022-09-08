# Impliments the Two Component Covariate Moderated EBNM


# impliment normal and point mass component distributions

.clamp <- function(v, .upper = 100, .lower = -100) {
  return(pmax(.lower, pmin(v, .upper)))
}

#' Convolved log density
#' @param  dist is a list with parameters for the normal component
#' @param betahat vector of oberveations
#' @param se vector of standard errors of observations
#' @return A vector of log densities log p(\hat\beta | se) = log N \hat\beta ; \mu, sigma^2 + se^2)
convolved_logpdf.normal <- function(dist, betahat, se) {
  sd <- sqrt(se^2 + dist$var)
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 100, -100)
  return(logp)
}

convolved_logpdf.point <- function(dist, betahat, se) {
  # return(dnorm(betahat, sd=se, log=T))
  sd <- se
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 1e3, -1e3)
  return(logp)
}

# data = list(bethat, se, X, Z)
# fit = list(data, logreg, f0, f1)

compute_post_assignment <- function(fit) {
  # ln p(z=1)/p(z=0)
  logit_pi <- compute_Xb.binsusie(fit)
  f0_loglik <- convolved_logpdf.point(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf.normal(fit$f1, fit$data$betahat, fit$data$se)
  logits <- f1_loglik - f0_loglik + logit_pi
  post_assignment <- sigmoid(logits)
  return(post_assignment)
}

compute_assignment_entropy <- function(fit) {
  p <- fit$data$y
  entropy <- -1 * (p * log(p) + (1 - p) * log(1 - p))
  return(entropy)
}

#' E[logp(y | f0, f1, z)]
#'
loglik.twococomo <- function(fit) {
  f0_loglik <- convolved_logpdf.point(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf.normal(fit$f1, fit$data$betahat, fit$data$se)
  loglik <- (1 - fit$data$y) * f0_loglik + fit$data$y * f1_loglik
  return(loglik)
}

compute_elbo.twococomo <- function(fit) {
  data_loglik <- sum(loglik.twococomo(fit))
  assignment_entropy <- sum(compute_assignment_entropy(fit), na.rm = TRUE)
  logreg_elbo <- compute_elbo.binsusie(fit)
  elbo <- data_loglik + assignment_entropy + logreg_elbo
  return(elbo)
}

init.twococomo <- function(data, L = 5, prior_mean = 0, prior_variance = 1, prior_weights = NULL) {
  # TODO:  check data has betahat, se, X, Z, N, etc
  data$X2 <- data$X^2

  # initialize component distribution
  f0 <- list(mu = 0)
  f1 <- list(mu = 0, var = 1)

  # initialize posterior assignment
  f0_loglik <- convolved_logpdf.point(f0, data$betahat, data$se)
  f1_loglik <- convolved_logpdf.normal(f1, data$betahat, data$se)
  logits <- f1_loglik - f0_loglik + 0 # TODO: set default prior odds to something other than 0

  data$y <- sigmoid(logits) # initialize component assignment
  fit <- init.binsusie(data,
    L = L,
    prior_mean = prior_mean,
    prior_variance = prior_variance,
    prior_weights = prior_weights
  )

  # add 2como specific attributes
  fit$f0 <- f0
  fit$f1 <- f1

  return(fit)
}

iter.twococomo <- function(fit, intercept = TRUE, estimate_prior_variance = TRUE) {
  # update assignments
  fit$data$y <- compute_post_assignment(fit)
  # update SuSiE
  fit <- iter.binsusie(fit, intercept, estimate_prior_variance)

  # TODO: update mixture components

  return(fit)
}

fit.twococomo <- function(data, maxiter = 100, tol = 1e-3) {
  fit <- init.twococomo(data)

  for (i in 1:maxiter) {
    fit <- iter.twococomo(fit)
    fit$elbo <- c(fit$elbo, compute_elbo.twococomo(fit))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
}


#' Point Normal SuSiE
#' Fits Point Normal SuSiE, observations are drawn from a two component mixture
#' (point-mass at zero, normal).
#' The prior log odds of being non-zero are a linear function of covariates X
#' and the effects have a SuSiE prior
#'  @param X a n x p matrix of covariates
#' @param betahat an n vector of observations
#' @param se a n vector of standard errors
#' @param Z fixed effect covaraites (including intercept). If null just a n x 1 matrix of ones
#' @param prior_mean the prior mean of each non-zero element of b. Either a scalar or vector of length L.
#' @param prior_variance the prior variance of each non-zero element of b. Either a scalar or vector of length L. If `estimate_prior_variance=TRUE` the value provides an initial estimate of prior variances
#' @param prior_weights prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix
#' @param intercept
#' @param estimate_prior_variance
#' @param s_init a logistic susie object to initialize with, NOTE if non-null, we ignore `prior_mean`, `prior_variance`, and `prior_weights`
#' @param returns a fit Binomial SuSiE model, which is compatable with summary functions from `susieR` package
#' @export
pointnormalsusie <- function(X,
                             betahat,
                             se,
                             Z = NULL,
                             L = min(10, ncol(X)),
                             prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
                             prior_variance = 1.0, # TODO: make scaled prior variance?
                             prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect... = hypers$pi
                             intercept = TRUE,
                             estimate_prior_variance = TRUE,
                             s_init = NULL, # previous fit with which to initialize NO
                             coverage = 0.95,
                             min_abs_corr = 0.5,
                             max_iter = 100,
                             tol = 0.001,
                             verbose = FALSE,
                             n_purity = 100) {
  ##################
  # checking / setup

  # Make data list
  data <- list(X = X, Z = Z, betahat = betahat, se = se)
  data$N <- 1 # binary observations
  # Initialize object
  if (is.null(s_init)) {
    fit <- init.twococomo(
      data,
      L = L,
      prior_mean = prior_mean,
      prior_variance = prior_variance,
      prior_weights = prior_weights
    )
  } else {
    fit <- s_init
  }

  # model fitting
  for (i in 1:max_iter) {
    fit <- iter.twococomo(fit, intercept, estimate_prior_variance)
    fit$elbo <- c(fit$elbo, compute_elbo.twococomo(fit))
    if (.converged(fit, tol)) {
      break
    }
  }

  fit <- binsusie_wrapup(fit)

  class(fit) <- c("binsusie", "susie")
  # 1. compute credible sets
  return(fit)
}
