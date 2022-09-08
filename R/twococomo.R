# Impliments the Two Component Covariate Moderated EBNM



convolved_logpdf.point <- function(dist, betahat, se) {
  # return(dnorm(betahat, sd=se, log=T))
  sd <- se
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 1e3, -1e3)
  return(logp)
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
