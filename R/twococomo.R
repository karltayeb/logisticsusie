# Impliments the Two Component Covariate Moderated EBNM


# impliment normal and point mass component distributions

.clamp <- function(v, .upper=100, .lower=-100){
  return(pmax(.lower, pmin(v, .upper)))
}

#' Convolved log density
#' @param  dist is a list with parameters for the normal component
#' @param betahat vector of oberveations
#' @param se vector of standard errors of observations
#' @return A vector of log densities log p(\hat\beta | se) = log N \hat\beta ; \mu, sigma^2 + se^2)
convolved_logpdf.normal <- function(dist, betahat, se){
  sd <- sqrt(se^2 + dist$var)
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log=TRUE)
  logp <- .clamp(logp, 100, -100)
  return(logp)
}

convolved_logpdf.point <- function(dist, betahat, se){
  # return(dnorm(betahat, sd=se, log=T))
  sd <- se
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log=TRUE)
  logp <- .clamp(logp, 1e3, -1e3)
  return(logp)
}

# data = list(bethat, se, X, Z)
# fit = list(data, logreg, f0, f1)

compute_post_assignment <- function(fit){
  # ln p(z=1)/p(z=0)
  logit_pi <- compute_Xb.binsusie(fit$logreg)
  f0_loglik <- convolved_logpdf.point(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf.normal(fit$f1, fit$data$betahat, fit$data$se)
  logits <- f1_loglik - f0_loglik + logit_pi
  post_assignment <- sigmoid(logits)
  return(post_assignment)
}

compute_assignment_entropy <- function(fit){
  p <- fit$post_assignment
  entropy = - 1 * (p * log(p) + (1-p) * log(1-p))
  return(entropy)
}

#' E[logp(y | f0, f1, z)]
#'
loglik.twococomo <- function(fit){
  f0_loglik <- convolved_logpdf.point(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf.normal(fit$f1, fit$data$betahat, fit$data$se)
  loglik <- (1 - fit$post_assignment) * f0_loglik + fit$post_assignment * f1_loglik
  return(loglik)
}

compute_elbo.twococomo <- function(fit){
  data_loglik <- sum(loglik.twococomo(fit))
  assignment_entropy <- sum(compute_assignment_entropy(fit), na.rm = TRUE)
  logreg_elbo <- compute_elbo.binsusie(fit$logreg)
  elbo <- data_loglik + assignment_entropy + logreg_elbo
  return(elbo)
}

init.twococomo <- function(data){
  #TODO:  check data has betahat, se, X, Z, N, etc

  # initialize component distribution
  f0 <- list(mu=0)
  f1 <- list(mu=0, var=1)

  # initialize posterior assignment
  f0_loglik <- convolved_logpdf.point(f0, data$betahat, data$se)
  f1_loglik <- convolved_logpdf.normal(f1, data$betahat, data$se)
  logits <- f1_loglik - f0_loglik + 0  # TODO: set default prior odds to something other than 0
  post_assignment <- sigmoid(logits)


  # initialize logistic regression
  logreg_data <- list(
    y = post_assignment,
    X = data$X,
    Z = data$Z,
    N = data$N
  )
  logreg <- fit.binsusie(logreg_data, maxiter = 2)

  fit <- list(
    data = data,
    f0=f0,
    f1=f1,
    post_assignment = post_assignment,
    logreg = logreg,
    elbo = c(-Inf)
  )
  return(fit)
}

iter.twococomo <- function(fit){
  # update mixture model
  fit$post_assignment <- compute_post_assignment(fit)
  fit$logreg$data$y <- fit$post_assignment
  fit$logreg <- iter.binsusie(fit$logreg)

  # TODO: update mixture components

  return(fit)
}

fit.twococomo <- function(data, maxiter=100, tol=1e-3){
  fit <- init.twococomo(data)

  for(i in 1:maxiter){
    fit <- iter.twococomo(fit)
    fit$elbo <- c(fit$elbo, compute_elbo.twococomo(fit))
    if (.converged(fit, tol)){
      break
    }
  }
  return(fit)
}
