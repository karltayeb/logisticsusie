# Generate X matrix ----
sim_X <- function(n = 1000, p = 50, length_scale = 50) {
  mix <- exp(-1 / length_scale * outer(seq(p), seq(p), "-")**2)
  X <- matrix(rnorm(n * p), nrow = n) %*% mix
  return(X)
}

sim_X_sparse <- function(n = 1000, p = 50, pi1 = 0.2, p_flip = 0.01) {
  X <- list()
  X[[1]] <- rbinom(n, 1, pi1)
  for (i in seq(2, p)) {
    x <- X[[i - 1]]
    on <- which(x == 1)
    off <- which(x == 0)

    flip <- rbinom(1, min(length(on), length(off)), p_flip)
    if (flip > 0) {
      x[sample(on, flip)] <- 0
      x[sample(off, flip)] <- 1
    }
    X[[i]] <- x
  }
  X <- Matrix::Matrix(do.call(cbind, X), sparse = T)
  return(X)
}


# Generate effects ----

sim_b_constant <- function(b = 1, p = 50, L = 3) {
  idx <- seq(L) * 10
  beta <- rep(b, L)
  return(list(beta = beta, idx = idx))
}

#' @export
sim_ser <- function(n = 1000,
                    p = 50,
                    N = 1,
                    idx = 1,
                    beta = 1,
                    beta0 = -2,
                    fsimX = sim_X,
                    fsimX_control = list()) {
  X <- rlang::exec(fsimX, n = n, p = p, !!!fsimX_control)
  Z <- matrix(rep(1, n), nrow = n)

  logits <- beta0 + beta * X[, idx]
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N, logits = logits
  )
  return(data)
}

#' simulate SER with 3 covariates
#' @export
sim_ser_with_covariates <- function(n = 1000, p = 50, N = 1, idx = 1) {
  X <- sim_X(n, p)
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
  X <- sim_X(n, p)
  Z <- matrix(rep(1, n), nrow = n)

  beta <- seq(L) * .2
  logits <- Matrix::drop(X[, seq(L) * 10] %*% beta)
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


sim_susie_sparse <- function(n = 1000, p = 50, L = 3, N = 1, pi1 = 0.2, transition = 0.8) {
  X <- sim_X_sparse(n, p, pi1, transition)
  Z <- matrix(rep(1, n), nrow = n)

  beta <- 3.
  logits <- -2 + beta * Matrix::rowSums(X[, 1:L])
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N, logits = logits
  )
  return(data)
}

sim_twococomo_sparse <- function(n = 1000, p = 50, L = 3, N = 1) {
  sim <- sim_susie_sparse(n, p, L, N)
  sim$beta <- 10 * rnorm(n) * sim$y
  sim$se <- 0.1 + rgamma(n, shape = 0.5)
  sim$betahat <- sim$beta + rnorm(n) * sim$se

  class(sim) <- "data_mococomo"
  return(sim)
}


sim_mn_susie <- function(n = 1000, p = 50, L = 3, N = 1, K = 10) {
  X <- sim_X(n = n, p = p)
  Z <- matrix(rep(1, n), nrow = n)

  Beta <- matrix(rnorm(p * K), nrow = p) / 10
  logits <- X %*% Beta

  Y <- t(do.call(cbind, purrr::map(seq(n), ~ rmultinom(1, 50, softmax(logits[.x, ])))))

  data <- list(
    X = X, Z = Z, y = Y, logits = logits
  )
  return(data)
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
