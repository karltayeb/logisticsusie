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

sim_y_ser <- function(X, beta0, beta, idx = NULL, N = 1) {
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (is.null(idx)) {
    idx <- sample(p, 1)
  }

  logits <- beta0 + beta * X[, idx]
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(y = y, logits = logits, N = N, beta = beta, beta0 = beta0, idx = idx)
  return(data)
}

sim_y_susie <- function(X, beta0, beta, re_var = 0, idx = NULL, N = 1) {
  n <- dim(X)[1]
  p <- dim(X)[2]

  if (length(beta) >= p) {
    stop("length(beta) must be less than number of columns of X")
  }
  if (is.null(idx)) {
    idx <- sample(p, length(beta))
  }

  logits <- beta0 + rowSums(beta * X[, idx, drop = F]) + (rnorm(n) * sqrt(re_var))
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)
  data <- list(y = y, logits = logits, N = N, beta = beta, beta0 = beta0, idx = idx)
  return(data)
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


sim_ser_with_random_effects <- function(n = 1000,
                                        p = 50,
                                        N = 1,
                                        idx = 1,
                                        beta = 1,
                                        beta0 = -2,
                                        re_var = 1,
                                        fsimX = sim_X,
                                        fsimX_control = list()) {
  X <- rlang::exec(fsimX, n = n, p = p, !!!fsimX_control)
  Z <- matrix(rep(1, n), nrow = n)

  logits <- beta0 + beta * X[, idx] + rnorm(n, sd = sqrt(re_var))
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X = X, Z = Z, y = y, N = N, logits = logits
  )
  return(data)
}


#' @title  Simulate data under the logistic SuSiE model
#' @description  Simulate data under the logistic SuSiE model
#' @param n numeric sample size
#' @param p numeric number of observed covariates
#' @param L numeric number of causal covariates (should be lower than p)
#' @param N numeric number of draw in the bionmial considered (set to 1 as default)
#' @param beta numeric if of length 1 assume that ervery effect have same beta coefficietns
#' otherwise should be  of length L, if missing beta <- seq(L) * .2
#' @param alpha numeric intercept in the lostic regression (control sparsity)
#' @export

sim_susie <- function(n = 1000, p = 50, L = 3, N = 1, beta, alpha = 0) {
  X <- sim_X(n, p)
  Z <- matrix(rep(1, n), nrow = n)
  if (missing(beta)) {
    beta <- seq(L) * .2
  }
  if (length(beta) == 1) {
    beta <- rep(beta, L)
  }


  logits <- alpha + Matrix::drop(X[, seq(L) * 10] %*% beta)
  p <- sigmoid(logits)
  y <- rbinom(n, N, p)

  data <- list(
    X         = X,
    Z         = Z,
    y         = y,
    N         = N,
    logits    = logits,
    effect    = beta,
    intercept = alpha
  )
  return(data)
}


#' @title  Simulate data under the twococomo SuSiE model
#' @description  Simulate data under the twococomo SuSiE model
#' @param n numeric sample size
#' @param p numeric number of observed covariates
#' @param L numeric number of causal covariates (should be lower than p)
#' @param N numeric number of draw in the bionmial considered (set to 1 as default)
#' @param beta numeric if of length 1 assume that ervery effect have same beta coefficietns
#' otherwise should be  of length L, if missing beta <- seq(L) * .2
#' @export
sim_twococomo <- function(n = 1000, p = 50, L = 3, N = 1, beta, alpha = 0) {
  if (missing(beta)) {
    sim <- sim_susie(n, p, L, N, alpha = alpha)
  } else {
    sim <- sim_susie(n, p, L, N, beta, alpha = alpha)
  }

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
sim_mococomo <- function(n = 1000, p = 50, L = 3, N = 1, beta, alpha = 0) {
  if (missing(beta)) {
    sim <- sim_susie(n, p, L, N, alpha = alpha)
  } else {
    sim <- sim_susie(n, p, L, N, beta, alpha = alpha)
  }
  sim$scales <- cumprod(c(1, rep(sqrt(2), 5)))
  sim$beta <- rnorm(n) * sim$y * sample(sim$scales, size = n, replace = T)
  sim$se <- 0.1 + rgamma(n, shape = 0.5)
  sim$betahat <- sim$beta + rnorm(n) * sim$se
  class(sim) <- "data_mococomo"
  return(sim)
}


#' @export
sim_ordinal_mod <- function(n = 1000,
                            p = 30,
                            p_act = 1,
                            se = 1,
                            n_class = 5,
                            beta_size = 1, # 1 moderatly informative , 1.5 strongly infrmative
                            alpha_start = 3.5, # 3.5, 2.5, 1.5 as start corresponds to sparse medium and dense signals
                            grid_s,
                            max_grid = 2, # represent signal power the larger the clearer signals
                            full_res = TRUE) {
  if (missing(grid_s)) {
    grid_s <- seq(0, max_grid, length.out = n_class)
  }
  if (missing(se)) {
    se <- runif(n)
  }
  if (length(se) == 1) {
    se <- rep(se, n)
  }


  beta <- sample(c(-beta_size, beta_size), size = p_act, replace = TRUE)



  alpha <- alpha_start + seq(0, 3, length.out = n_class)
  # res_summary= length(which(cebmn_fdr<0.05))/n
  # length(which(camt_lfdr<0.05))/n
  # length(which(ihw_fdr<0.05))/n

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  # simulation under ordinal model
  if (length(beta) == 1) {
    pred <- sigmoid(matrix(alpha,
      byrow = TRUE,
      nrow = n,
      ncol = n_class
    ) + matrix(rep(X[, c(1:p_act)] * (beta), n_class),
      ncol = n_class
    ))
  } else {
    pred <- sigmoid(matrix(alpha,
      byrow = TRUE,
      nrow = n,
      ncol = n_class
    ) + matrix(rep(X[, c(1:p_act)] %*% (beta), n_class),
      ncol = n_class
    ))
  }



  ind_prob <- matrix(NA, ncol = n_class, nrow = n)
  for (i in 1:nrow(ind_prob)) {
    for (j in 1:ncol(ind_prob)) {
      if (j == 1) {
        ind_prob[i, j] <- pred[i, j]
      } else {
        ind_prob[i, j] <- pred[i, j] - pred[i, (j - 1)]
      }
    }
  }
  ind_prob <- ind_prob / apply(ind_prob, 1, sum)

  obs <- rep(NA, n)
  true_obs <- rep(NA, n)
  class <- rep(NA, n)
  for (i in 1:n)
  {
    class[i] <- sample(size = 1, 1:n_class, prob = ind_prob[i, ])
    true_obs[i] <- ifelse(class[i] == 1, 0, rnorm(1, 0, grid_s[class[i]]))
    obs[i] <- true_obs[i] + rnorm(1, sd = se[i])
  }
  out <- list(
    true_obs = true_obs,
    obs = obs,
    X = X,
    se = se,
    class = class,
    sim_param = c(n, p, p_act, se, n_class)
  )
  return(out)
}
