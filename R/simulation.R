#' @export

sim_ser <- function(n = 1000, p = 50, N = 1, idx = 1) {
  mix <- exp(-0.1 * outer(seq(p), seq(p), "-")**2)
  X <- matrix(rnorm(n * p), nrow = n) %*% mix
  Z <- matrix(rep(1, n), nrow = n)

  beta <- 1.
  p <- sigmoid(X[, idx])
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N
  )
  return(data)
}

#' simulate SER with 3 covariates
#' @export
sim_ser_with_covariates <- function(n = 1000, p = 50, N = 1, idx = 1) {
  mix <- exp(-0.1 * outer(seq(p), seq(p), "-")**2)
  X <- matrix(rnorm(n * p), nrow = n) %*% mix
  Z <- matrix(rnorm(n * 3), nrow = n)

  p <- sigmoid(X[, idx] + rowSums(Z))
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N
  )
  return(data)
}
#' @export

sim_susie <- function(n = 1000, p = 50, L = 3, N = 1) {
  mix <- exp(-0.1 * outer(seq(p), seq(p), "-")**2)
  X <- matrix(rnorm(n * p), nrow = n) %*% mix # mix so that there are correlated features
  X <- X[, sample(p)] # permute so that adjacent columns aren't very correlated
  Z <- matrix(rep(1, n), nrow = n)

  beta <- 1.
  logits <- rowSums(X[, 1:L])
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N, logits = logits
  )
  return(data)
}

sim_X_sparse <- function(n = 1000, p = 50, pi1 = 0.2, transition = 0.8) {
  X <- list()
  X[[1]] <- rbinom(n, 1, pi1)
  for (i in seq(2, p)) {
    x <- X[[i - 1]]
    X[[i]] <- ifelse(rbinom(n, 1, transition) == 1, x, abs(x - 1))
  }
  X <- Matrix::Matrix(do.call(cbind, X), sparse = T)
  return(X)
}

sim_susie_sparse <- function(n = 1000, p = 50, L = 3, N = 1, pi1 = 0.2, transition = 0.8) {
  X <- sim_X_sparse(n, p, pi1, transition)
  Z <- matrix(rep(1, n), nrow = n)

  beta <- 1.
  logits <- Matrix::rowSums(X[, 1:L])
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N, logits = logits
  )
  return(data)
}

#' @export
sim_twococomo <- function(n = 1000, p = 50, L = 3, N = 1) {
  sim <- sim_susie(n, p, L, N)
  sim$beta <- rnorm(n) * sim$y
  sim$se <- 0.1 + rgamma(n, shape = 0.5)
  sim$betahat <- sim$beta + rnorm(n) * sim$se

  class(sim) <- "data_mococomo"
  return(sim)
}

#' @export
sim_mococomo <- function(n = 1000, p = 50, L = 3, N = 1) {
  sim <- sim_susie(n, p, L, N)
  sim$scales <- cumprod(c(1, rep(sqrt(2), 5)))
  sim$beta <- rnorm(n) * sim$y * sample(sim$scales, size = n, replace = T)
  sim$se <- 0.1 + rgamma(n, shape = 0.5)
  sim$betahat <- sim$beta + rnorm(n) * sim$se
  class(sim) <- "data_mococomo"
  return(sim)
}
