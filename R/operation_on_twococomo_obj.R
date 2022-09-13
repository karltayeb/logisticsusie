
# data = list(bethat, se, X, Z)
# fit = list(data, logreg, f0, f1)

compute_post_assignment <- function(fit) {
  # ln p(z=1)/p(z=0)
  logit_pi <- compute_Xb.binsusie(fit)
  f0_loglik <- convolved_logpdf(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf(fit$f1, fit$data$betahat, fit$data$se)
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
  f0_loglik <- convolved_logpdf(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf(fit$f1, fit$data$betahat, fit$data$se)
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

#' Initialize Two CoCoMo
#' @param data list of data for twococomo
#' @param L number of single effects
#' @param f1_dist the name of `component_distribution` for the non-null component
#' @param f1_params list of parameters to initialize f1, must be a subset of those accepted by the constructor for f1
#' @param prior_mean the prior mean of each non-zero element of b. Either a scalar or vector of length L.
#' @param prior_variance the prior variance of each non-zero element of b. Either a scalar or vector of length L. If `estimate_prior_variance=TRUE` the value provides an initial estimate of prior variances
#' @param prior_weights prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix
#' @param returns am initialized twococomo object
init.twococomo <- function(data,
                           L = 5,
                           f1_dist = "normal",
                           f1_params = list(),
                           prior_mean = 0,
                           prior_variance = 1,
                           prior_weights = NULL) {
  # TODO:  check data has betahat, se, X, Z, N, etc
  data$X2 <- data$X^2

  # initialize component distribution
  f0 <- point_component(mu = 0)

  f1 <- rlang::exec(
    paste0(f1_dist, "_component"), # constructor name
    !!!f1_params # unpack list of params
  )

  # initialize posterior assignment
  f0_loglik <- convolved_logpdf(f0, data$betahat, data$se)
  f1_loglik <- convolved_logpdf(f1, data$betahat, data$se)
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

iter.twococomo <- function(fit, intercept = TRUE, estimate_prior_variance = TRUE, estimate_f1 = TRUE) {
  # update assignments
  fit$data$y <- compute_post_assignment(fit)

  # update SuSiE
  fit <- iter.binsusie(fit, intercept, estimate_prior_variance)

  # TODO: update mixture components
  if (estimate_f1) {
    fit$f1 <- update_params(
      fit$f1, fit$data$betahat, fit$data$se,
      weights = fit$data$y
    )
  }
  return(fit)
}

fit.twococomo <- function(data, maxiter = 100, tol = 1e-3, estimate_f1 = TRUE) {
  fit <- init.twococomo(data)

  for (i in 1:maxiter) {
    fit <- iter.twococomo(fit, estimate_f1 = estimate_f1)
    fit$elbo <- c(fit$elbo, compute_elbo.twococomo(fit))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
}
