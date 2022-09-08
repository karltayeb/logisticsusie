
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
