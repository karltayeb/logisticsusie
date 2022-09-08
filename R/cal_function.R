
#' @param gamma cumulative probability p(z <= k)
#' @param pi_tilde conditional probability p(z = k+1 | z > k)
#' @return cumulative probability p(z <= k + 1)
.accumulate_gamma <- function(gamma, pi_tilde) {
  pi <- pi_tilde * (1. - gamma)
  gamma <- gamma + pi
  return(gamma)
}

#' Convert \tilde pi to \pi
#' Converts conditional probabilities p(z = k | z >=k) k = 1... K
#' to marginal probilities p(z =k) k = 1... K
#' @param pi_tilde vector of conditional probabilities (except for last one p(z=K | Z >=K) = 1)
#' @return pi a vector of probabilities summing to one
.pi_tilde2pi <- function(pi_tilde) {
  gamma <- purrr::accumulate(pi_tilde, .accumulate_gamma)
  pi <- c(gamma[1], tail(gamma, -1) - head(gamma, -1), 1 - tail(gamma, 1))
  return(pi)
}


#' Predict to log pi
#' Computes E[log p(z | \beta, \omega)] from E[Xb_1], ... E[Xb_{K-1}]
#' Where expectations are taken under variational approximation
#' @param xb vector, length K-1 of conditional log odds for first K-1 components
#' @return a vector of length K of (expected) log prior assignment probabilities
#' @example
#' exp(.predict2logpi(rep(0, 3)))  # = c(0.5, 0.25, 0.25)
.predict2logpi <- function(xb) {
  K <- length(xb) + 1
  xbcum <- cumsum(xb)
  kln2 <- seq(K - 1) * log(2)
  logpi <- -kln2 + xb - 0.5 * xbcum # un-normalized k=1, ..., K-1
  logpiK <- -tail(kln2, 1) - 0.5 * tail(xbcum, 1) # for component K
  logpi <- c(logpi, logpiK)
  logpi <- logpi - logSumExp(logpi)
  return(logpi)
}

#' Compute Log Prior Assignment Probabilities
#' For each data point return covariate-dependent prior mixture probabilities
#' @param fit a MoCoCoMo fit object
#' @return an n x K matrix of log prior probabilities for each data point
.compute_prior_assignment <- function(xb) {
  .predict2logpi(xb)
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







#' Compute Data Log Likelihood
#' Compute p(betahat_i| z=k, se_i) i = 1..., n, k=1... K
#' @param fit an MoCoCoMo fit object
#' @return n x K matrix of  log likelihood of each data point for each component
compute_data_loglikelihood <- function(fit) {
  # compute loglikelihood of each data point under each component distribution
  # TODO: replace with generic convolved_logpdf function
  data_loglik <- do.call(cbind, purrr::map(
    fit$f_list, ~ convolved_logpdf.normal(.x, fit$data$betahat, fit$data$se)
  ))
  return(data_loglik)
}


.asignment_loglik <- function(xb, xb2) {
  K <- length(xb) + 1
  xbcum <- cumsum(xb)
  kln2 <- seq(K - 1) * log(2)
  logpi <- -kln2 + xb - 0.5 * xbcum # un-normalized k=1, ..., K-1
  logpiK <- -tail(kln2, 1) - 0.5 * tail(xbcum, 1) # for component K
  logpi <- c(logpi, logpiK)
  logpi <- logpi - logSumExp(logpi)
  return(logpi)
}



#' polya Gamma Mean
#' Compute the mean of a Polya Gamma RV w ~ PG(b, c)
#'  @param b a PG parameter
#'  @param c a PG parameter
#'  @return E[w] the mean of w ~ PG(b, c)
pg_mean <- function(b, c) {
  mu <- 0.5 * b / c * tanh(c / 2)

  # deal with case of c = 0 mean is b/4
  idx <- is.na(mu)
  mu[idx] <- b[idx] / 4
  return(mu)
}

####
# KLs
###

#' Normal KL
#' Compute KL(N(mu, var) || N(mu0, var0))
#' @param mu first mean parameter of
#' @param var fist variance parameter
#' @param mu0 second mean parameter
#' @param var0 second variance parameter
#' @return KL(N(mu, var) || N(mu0, var0))
normal_kl <- function(mu, var, mu0 = 0., var0 = 1.) {
  kl <- 0.5 * (log(var0) - log(var) + var / var0 + (mu - mu0)^2 / var0 - 1)
  return(kl)
}

#' Polya-Gamma KL
#' Compute KL(PG(b, c) || N(b, 0))
#' @param b first mean parameter of
#' @param c fist variance parameter
#' @return Return KL-divergence KL(PG(b, c) || N(b, 0))
pg_kl <- function(b, c) {
  kl <- -0.5 * c^2 * pg_mean(b, c) + b * log(cosh(c / 2))
  return(kl)
}

#' Categorical KL
#' Compute KL(Categorical(alpha) || Categorical(pi))
#' @param alpha vector of probabilities that sum to 1
#' @param pi vector of probabilities that sum to 1
#' @return Return KL-divergence KL(Categorical(alpha) || Categorical(pi))
categorical_kl <- function(alpha, pi) {
  idx <- alpha > 1 / (length(alpha) * 1e3)
  sum((alpha * (log(alpha) - log(pi))), na.rm = TRUE)
}


#' Categorical Entropy
#' Compute entropy of Categorical(pi) = -E[log p]
#' @param pi vector of probabilities that sum to 1
#' @return Return entropy -Elog(x)
categorical_entropy <- function(pi) {
  idx <- pi > 1 / (length(pi) * 1e3)
  entropy <- -1 * sum(pi * log(pi), na.rm = TRUE)
}
