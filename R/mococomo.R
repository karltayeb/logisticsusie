# Impliments covariate moderated ASH "MOre COmponents COvariate MOderated"

#' @param gamma cumulative probability p(z <= k)
#' @param pi_tilde conditional probability p(z = k+1 | z > k)
#' @return cumulative probability p(z <= k + 1)
.accumulate_gamma <- function(gamma, pi_tilde){
  pi = pi_tilde * (1. - gamma)
  gamma = gamma + pi
  return(gamma)
}

#' Convert \tilde pi to \pi
#' Converts conditional probabilities p(z = k | z >=k) k = 1... K
#' to marginal probilities p(z =k) k = 1... K
#' @param pi_tilde vector of conditional probabilities (except for last one p(z=K | Z >=K) = 1)
#' @return pi a vector of probabilities summing to one
.pi_tilde2pi <- function(pi_tilde){
  gamma <- purrr::accumulate(pi_tilde, .accumulate_gamma)
  pi <-  c(gamma[1], tail(gamma, -1) - head(gamma, -1), 1 - tail(gamma, 1))
  return(pi)
}


#' Predict to log pi
#' Computes E[log p(z | \beta, \omega)] from E[Xb_1], ... E[Xb_{K-1}]
#' Where expectations are taken under variational approximation
#' @param xb vector, length K-1 of conditional log odds for first K-1 components
#' @return a vector of length K of (expected) log prior assignment probabilities
#' @example
#' exp(.predict2logpi(rep(0, 3)))  # = c(0.5, 0.25, 0.25)
.predict2logpi <- function(xb){
  K <- length(xb) + 1
  xbcum <- cumsum(xb)
  kln2 <- seq(K - 1) * log(2)
  logpi <- -kln2 + xb - 0.5 * xbcum  # un-normalized k=1, ..., K-1
  logpiK <- -tail(kln2, 1) - 0.5 * tail(xbcum, 1)  # for component K
  logpi <- c(logpi, logpiK)
  logpi <- logpi - logSumExp(logpi)
  return(logpi)
}

#' Compute Log Prior Assignment Probabilities
#' For each data point return covariate-dependent prior mixture probabilities
#' @param fit a MoCoCoMo fit object
#' @return an n x K matrix of log prior probabilities for each data point
.compute_prior_assignment <- function(xb){
  .predict2logpi(xb)
}


#' Get N x K-1 matrix of PG variational parameters
#' TODO: deprecate this when we call update functions directly rather than
#' passing to logistic susie subroutine
get_xi.mococomo <- function(fit){
  xi <- do.call(cbind, purrr::map(fit$logreg_list, ~purrr::pluck(.x, 'params', 'xi')))  # N x K-1
  return(xi)
}

compute_Xb.mococomo <- function(fit){
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, compute_Xb.binsusie))  # N x K-1
  return(Xb)
}

compute_prior_assignment <- function(fit){
  # TODO alias compute_Xb with predict so that it works with other functions?
  # TODO make sure GLM predict outputs log-odds scale?
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, compute_Xb.binsusie))  # N x K-1
  res <- do.call(rbind, apply(Xb, 1, .predict2logpi, simplify = F))  # N x K
  return(res)
}


.compute_prior_assignment2 <- function(xb){
  log(.pi_tilde2pi(sigmoid(xb)))
}

#' DONT USE, this ignores the PG agumentation and treats E[Xb_k] as the conditional log odds of sampling component k
#' and simply converts these quantitities to a vector of probabilities
compute_prior_assignment2 <- function(fit, log=TRUE){
  # TODO alias compute_Xb with predict so that it works with other functions?
  # TODO make sure GLM predict outputs log-odds scale?
  pi_tilde <- do.call(cbind, purrr::map(fit$logreg_list, ~sigmoid(compute_Xb.binsusie(.x))))  # N x K-1
  res <- do.call(rbind, apply(pi_tilde, 1, .pi_tilde2pi, simplify = F))  # N x K

  if(log){res <- log(res)}  # exponentiate if log=FALSE
  return(res)
}


#' Compute Data Log Likelihood
#' Compute p(betahat_i| z=k, se_i) i = 1..., n, k=1... K
#' @param fit an MoCoCoMo fit object
#' @return n x K matrix of  log likelihood of each data point for each component
compute_data_loglikelihood <- function(fit){
  # compute loglikelihood of each data point under each component distribution
  # TODO: replace with generic convolved_logpdf function
  data_loglik <- do.call(cbind, purrr::map(
    fit$f_list, ~convolved_logpdf.normal(.x, fit$data$betahat, fit$data$se)))
  return(data_loglik)
}


.asignment_loglik <- function(xb, xb2){
  K <- length(xb) + 1
  xbcum <- cumsum(xb)
  kln2 <- seq(K - 1) * log(2)
  logpi <- -kln2 + xb - 0.5 * xbcum  # un-normalized k=1, ..., K-1
  logpiK <- -tail(kln2, 1) - 0.5 * tail(xbcum, 1)  # for component K
  logpi <- c(logpi, logpiK)
  logpi <- logpi - logSumExp(logpi)
  return(logpi)
}

#' Compute Expected Assignment Log Likelihood
#' @return N x K matrix each row has E[p(z=k | \beta, \omega)] for k = 1,..., K
compute_assignment_loglikelihood.mococomo <- function(fit, normalize=TRUE){
  K <- length(fit$f_list)
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, compute_Xb.binsusie))  # N x K-1
  Xbcum <- do.call(rbind, apply(Xb, 1, cumsum, simplify = F))
  kln2 <- seq(K-1) * log(2)

  # second moments only appear in normalizing constant
  C <- 0
  if(normalize){
    Xb2 <- do.call(cbind, purrr::map(fit$logreg_list, compute_Xb2.binsusie))  # N x K - 1
    omega <- do.call(cbind, purrr::map(fit$logreg_list, compute_omega))  # N x K-1
    C <- -0.5 * rowSums(Xb2 * omega)
  }

  logp <- Xb - 0.5 * Xbcum  + C
  logp <- do.call(rbind, apply(logp, 1, function(x) x - kln2, simplify = F))

  logpK <- -(K-1) * log(2) - 0.5 * Xbcum[, K-1] + C
  logp <- cbind(logp, logpK)
  return(logp)
}



compute_assignment_jj_bound.mococomo <- function(fit){
  K <- length(fit$f_list)
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, compute_Xb.binsusie))  # N x K-1
  Xi <- do.call(cbind, purrr::map(fit$logreg_list, ~purrr::pluck(.x, 'params', 'xi')))  # N x K-1

  f <- function(xi, xb){
    tmp <- cumsum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb) + xb
    tmpK <- sum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb)
    jj <- c(tmp, tmpK)
    return(jj)
  }
  jj <- do.call(rbind, purrr::map(seq(nrow(Xb)), ~ f(Xi[.x,], Xb[.x,])))

  # JJ bound, should be faster with matrix operations?
  # tmp <- rowCumSum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * Xb)
  # jjK <- tmp[, K-1]
  # jj <- cbind(Xb + tmp, jjK)
  return(jj)
}

#' Compute Posterior Assignment Probabilities
#' For each data point return posterior assignment probabilities
#' @param fit a MoCoCoMo fit object
#' @return an n x K matrix of log posterior probabilities for each data point
compute_posterior_assignment <- function(fit, log=TRUE){
  data_loglik <- fit$data_loglik
  assignment_loglik <- compute_assignment_jj_bound.mococomo(fit)
  # assignment_loglik <- compute_assignment_loglikelihood.mococomo(fit)

  # normalize
  res <- do.call(
    rbind,
    apply(data_loglik + assignment_loglik, 1, function(x) x - logSumExp(x), simplify = F)
  )
  if(!log){res <- exp(res)}  # exponentiate if log=FALSE
  return(res)
}


compute_elbo3.mococomo <- function(fit){
  post_assignment <- fit$post_assignment
  data_loglik <- fit$data_loglik
  jj <- compute_assignment_jj_bound.mococomo(fit)

  kl_susie <- sum(purrr::map_dbl(fit$logreg_list, compute_kl.binsusie))
  assignment_entropy <- sum(apply(post_assignment, 1, categorical_entropy))

  elbo <- sum(post_assignment * (jj + data_loglik)) + assignment_entropy - kl_susie
  return(elbo)
}

compute_elbo2.mococomo <- function(fit){
  post_assignment <- fit$post_assignment
  assignment_loglik <- compute_assignment_loglikelihood.mococomo(fit)
  ll <- sum(post_assignment * (fit$data_loglik + assignment_loglik))
  assignment_entropy <- sum(apply(post_assignment, 1, categorical_entropy))

  kl_susie <- sum(purrr::map_dbl(fit$logreg_list, compute_kl.binsusie))
  kl_omega <- sum(purrr::map_dbl(fit$logreg_list, function(x) sum(pg_kl(x$data$N, x$params$xi))))
  kl <- kl_susie + kl_omega

  return (ll + assignment_entropy - kl)
}

compute_elbo.mococomo <- function(fit){
  K <- length(fit$f_list)
  N <- .expected_trials(fit$post_assignment)
  # MAKE SURE susie has access to the right posteriors!
  # and update omega so jj bound is tight
  for(k in seq(K-1)){
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- N[, k]
    fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit$logreg_list[[k]])
    fit$logreg_list[[k]]$params$tau <- compute_tau(fit$logreg_list[[k]])
  }


  # TODO: do we want to just store post_assignment after each update?
  post_assignment <- fit$post_assignment
  data_loglik <- fit$data_loglik

  ll <- sum(post_assignment * data_loglik)
  assignment_entropy <- sum(apply(post_assignment, 1, categorical_entropy))
  logreg_elbo <- sum(purrr::map_dbl(fit$logreg_list, compute_elbo.binsusie))
  elbo <- ll + assignment_entropy + logreg_elbo
  return(elbo)
}

#' Expected Binomial Trials
#' Compute expected  number of binomial trials
#' MoCoCoMo fits K-1 binomial regressions, this is the E-step for N
#' @param Y an N x K matrix of assignment probabilities
.expected_trials <- function(Y) {
  cumY <- do.call(rbind, apply(Y, 1, cumsum, simplify = F))
  N <- 1. - cumY + Y
}

init.mococomo <- function(data){
  # TODO check input data has X, Z, betahat, se, etc.
  # scale mixture of normals w/ zero mean-- TODO: make this initializable
  scales <- cumprod(c(1., rep(sqrt(2), 5)))
  f_list <- purrr::map(scales, ~list(mu=0, var=.x^2))

  K = length(scales)
  p <- ncol(data$X)
  n <- nrow(data$X)

  # initialize posterior assignment
  Y <- matrix(rep(1, K*n)/K, nrow=n)
  N <- .expected_trials(Y)

  data$y <- Y
  data$N <- N

  # initialize K-1 logistic SuSiE
  .init_logreg <- function(y, N){
    dat <- list(
      X = data$X,
      Z = data$Z,
      y = y,
      N = N
    )
    fit <- init.binsusie(dat)
    return(fit)
  }
  logreg_list <- purrr::map(seq(K-1), ~ .init_logreg(Y[, .x], N[, .x]))
  fit <- list(
    data = data,
    logreg_list=logreg_list,
    f_list=f_list,
    post_assignment=Y,
    N=N
  )

  # only need to do this once when component probabilities are fixed
  fit$data_loglik <- compute_data_loglikelihood(fit)
  return(fit)
}


iter.mococomo <- function(fit, update_assignment=T, update_logreg=T){
  K <- length(fit$f_list)

  # updates posterior assignments
  if(update_assignment){
    fit$post_assignment <- compute_posterior_assignment(fit, log=F)
    fit$N  <- .expected_trials(fit$post_assignment)
  }

  if(update_logreg){
    logreg_list <- list()
    for(k in seq(K-1)){
      logreg <- fit$logreg_list[[k]]

      # pass assignments to logreg
      logreg$data$y <- fit$post_assignment[, k]
      logreg$data$N <- fit$N[, k]

      # update logreg
      logreg_list[[k]] <- iter.binsusie(logreg)
    }
    fit$logreg_list <- logreg_list
  }

  return(fit)
}


is.odd <- function(x){x%%2 == 1}
is.even <- function(x){x%%2 == 0}

fit.mococomo <- function(data, maxiter=100, tol=1e-3){
  fit <- init.mococomo(data)
  fit$elbo <- compute_elbo.mococomo(fit)
  for(i in 1:maxiter){
    fit <- iter.mococomo(fit, is.even(i), is.odd(i))
    fit$elbo <- c(fit$elbo, compute_elbo3.mococomo(fit))

    # print(paste('asgn:', is.even(i), 'logreg:', is.odd(i), 'elbo: ', tail(fit$elbo, 1)))
    if (.converged(fit, tol)){
      break
    }
  }
  return(fit)
}

#
# data <- sim_twococomo()
# fit <- fit.mococomo(data, maxiter=20)
# plot(fit$elbo)
# .monotone(fit$elbo)






