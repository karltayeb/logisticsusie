
#' @param gamma cumulative probability p(z <= k)
#' @param pi_tilde conditional probability p(z = k+1 | z > k)
#' @return cumulative probability p(z <= k + 1)
accumulate_gamma <- function(gamma, pi_tilde) {
  pi <- pi_tilde * (1. - gamma)
  gamma <- gamma + pi
  return(gamma)
}

#' Convert \tilde pi to \pi
#'
#' Converts conditional probabilities p(z = k | z >=k) k = 1... K
#' to marginal probilities p(z =k) k = 1... K
#' @param pi_tilde vector of conditional probabilities (except for last one p(z=K | Z >=K) = 1)
#' @return pi a vector of probabilities summing to one
pi_tilde2pi <- function(pi_tilde) {
  gamma <- purrr::accumulate(pi_tilde, accumulate_gamma)
  pi <- c(gamma[1], tail(gamma, -1) - head(gamma, -1), 1 - tail(gamma, 1))
  return(pi)
}


#' Predict to log pi
#' Computes E[log p(z | \beta, \omega)] from E[Xb_1], ... E[Xb_{K-1}]
#' Where expectations are taken under variational approximation
#' @param xb vector, length K-1 of conditional log odds for first K-1 components
#' @return a vector of length K of (expected) log prior assignment probabilities
#' @example
#' exp(predict2logpi(rep(0, 3)))  # = c(0.5, 0.25, 0.25, 0.125)
predict2logpi <- function(xb) {
  K <- length(xb) + 1
  xbcum <- cumsum(xb)
  kln2 <- seq(K - 1) * log(2)
  logpi <- -kln2 + xb - 0.5 * xbcum # un-normalized k=1, ..., K-1
  logpiK <- -tail(kln2, 1) - 0.5 * tail(xbcum, 1) # for component K
  logpi <- c(logpi, logpiK)
  logpi <- logpi - logSumExp(logpi)
  return(logpi)
}

