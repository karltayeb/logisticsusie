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


sim_susie <- function(n=1000, p=50, L=3, N=1){
  mix <- exp(- 0.1 * outer(seq(p), seq(p), '-')**2)
  X <- matrix(rnorm(n*p), nrow=n) %*% mix  # mix so that there are correlated features
  X <- X[, sample(p)]  # permute so that adjacent columns aren't very correlated
  Z <- matrix(rep(1, n), nrow=n)

  beta <- 1.
  logits <- rowSums(X[, 1:L])
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z=Z, y=y, N=N, logits=logits
  )
  return(data)
}


sim_twococomo <- function(n=1000, p=50, L=3, N=1){
  sim <- sim_susie(n, p, L, N)
  sim$beta <- rnorm(n) * sim$y
  sim$se <- 0.1 + rgamma(n, shape=0.5)
  sim$betahat <- sim$beta + rnorm(n) * sim$se
  return(sim)
}

