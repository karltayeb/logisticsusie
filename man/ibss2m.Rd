% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_random_effect_vb.R
\name{ibss2m}
\alias{ibss2m}
\title{2 moments Iterative Bayesian Stepwise Selectioin}
\usage{
ibss2m(
  X,
  y,
  L = 10,
  prior_variance = 1,
  prior_weights = NULL,
  tol = 0.001,
  maxit = 100,
  estimate_intercept = T,
  estimate_prior_variance = F,
  track_elbo = F
)
}
\description{
2 moments Iterative Bayesian Stepwise Selectioin
}
