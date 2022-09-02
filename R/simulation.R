sim_ser <- function(n=1000, p=50, N=1, idx=1){
  mix <- exp(- 0.1 * outer(seq(p), seq(p), '-')**2)
  X <- matrix(rnorm(n*p), nrow=n) %*% mix
  Z <- matrix(rep(1, n), nrow=n)

  beta <- 1.
  p <- sigmoid(X[, idx])
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z=Z, y=y, N=N
  )
  return(data)
}

#' simulate SER with 3 covariates
sim_ser_with_covariates <- function(n=1000, p=50, N=1, idx=1){
  mix <- exp(- 0.1 * outer(seq(p), seq(p), '-')**2)
  X <- matrix(rnorm(n*p), nrow=n) %*% mix
  Z <- matrix(rnorm(n*3), nrow=n)

  p <- sigmoid(X[, idx] + rowSums(Z))
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z=Z, y=y, N=N
  )
  return(data)
}
