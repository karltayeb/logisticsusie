#' polya Gamma Mean
#' Compute the mean of a Polya Gamma RV w ~ PG(b, c)
#'  @param b a PG parameter
#'  @param c a PG parameter
#'  @return E[w] the mean of w ~ PG(b, c)
pg_mean <- function(b, c) {
  mu <- Matrix::drop(0.5 * b / c * tanh(c / 2))

  # deal with case of c = 0 mean is b/4
  idx <- is.na(mu) # does indexing work form matrix entries?
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
