
###
# Single Effect Regression
###

###
# L-indexed getters
###

#' idx is for indexing when we run susie
.get_alpha <- function(fit, idx=NULL){
  if(is.null(idx)){
    alpha <- fit$params$alpha
  } else{
    alpha <- fit$params$alpha[idx,]
  }
  return(alpha)
}

.get_mu <- function(fit, idx=NULL){
  if(is.null(idx)){
    mu <- fit$params$mu
  } else{
    mu <- fit$params$mu[idx,]
  }
  return(mu)
}

.get_var <- function(fit, idx=NULL){
  if(is.null(idx)){
    var <- fit$params$var
  } else{
    var <- fit$params$var[idx,]
  }
  return(var)
}

.get_pi <- function(fit, idx=NULL){
  if(is.null(idx)){
    pi <- fit$hypers$pi
  } else{
    pi <- fit$hypers$pi[idx,]
  }
  return(pi)
}

.get_delta <- function(fit, idx=NULL){
  if(is.null(idx)){
    delta <- fit$params$delta
  } else{
    delta <- fit$params$delta[idx,]
  }
  return(delta)
}

.get_mu0 <- function(fit, idx=NULL){
  if(is.null(idx)){
    mu0 <- fit$hypers$mu0
  } else{
    mu0 <- fit$hypers$mu0[idx]
  }
  return(mu0)
}

.get_var0 <- function(fit, idx=NULL){
  if(is.null(idx)){
    var0 <- fit$hypers$sigma0^2
  } else{
    var0 <- fit$hypers$sigma0[idx]^2
  }
  return(var0)
}

###
# K-indexed getters
###

.get_y <- function(fit, kidx=NULL){
  if(is.null(kidx)){
    y <- fit$data$y
  } else{
    y <- fit$data$y[, kidx]
  }
  return(y)
}

.get_N <- function(fit, kidx=NULL){
  if(is.null(kidx)){
    N <- fit$data$N
  } else{
    N <- fit$data$N[, kidx]
  }
  return(N)
}

.get_xi <- function(fit, kidx=NULL){
  if(is.null(kidx)){
    xi <- fit$params$xi
  } else{
    xi <- fit$params$xi[, kidx]
  }
  return(xi)
}


#' Ebeta.binser
#' Compute expected value of coefficients beta under approximate posterior
#' @param fit Binomial SER object
#' @return Return E[\beta]
Ebeta.binser <- function(fit, idx=NULL){
  b <- drop(.get_alpha(fit, idx) * .get_mu(fit, idx))
  return(b)
}

#' KL(q(b, \omega) || p(\beta, \omega))
.compute_ser_kl <- function(alpha, pi, mu, var, mu0, var0){
  kl <- categorical_kl(alpha, pi)
  kl <- kl + sum(alpha * normal_kl(
    mu, var, mu0, var0))
  return(kl)
}

#' compute_kl.binser
#' Compute KL divergence between approximate posterior and prior
#' @param fit Binomial SER object
#' @return Return KL(q(\beta) || p(\beta))
compute_kl.binser <- function(fit){
  return(.compute_ser_kl(
    .get_alpha(fit),
    .get_pi(fit),
    .get_mu(fit),
    .get_var(fit),
    .get_mu0(fit),
    .get_var0(fit)
  ))
}

#' fixed effects, intercept, shift
compute_Zd.binser <- function(fit, idx=NULL, shift=0){
  Zd <- drop(fit$data$Z %*% .get_delta(fit, idx)) + shift
  return(Zd)
}

#' compute_Xb.binser
#' Compute expected predicted log odds
#' @param fit Binomial SER object
#' @return Return E[Xb]
compute_Xb.binser <- function(fit, idx=NULL, shift=0){
  Xb <- drop(fit$data$X %*% Ebeta.binser(fit, idx))
  Zd <- compute_Zd.binser(fit, idx, shift)
  Xb <- Xb + Zd
  return(Xb)
}

#' compute_Xb2.binser
#' Compute second moments of predicted log odds
#' @param fit Binomial SER object
#' @return Return E[Xb]
compute_Xb2.binser <- function(fit, idx=NULL, shift=0){
  Xb <- drop(fit$data$X %*% Ebeta.binser(fit, idx))
  Zd <- compute_Zd.binser(fit, idx, shift)

  b2 <- .get_alpha(fit, idx) * (.get_mu(fit, idx)^2 + .get_var(fit, idx))
  Xb2 <- drop(fit$data$X^2 %*% b2)

  Xb2 <- Xb2 + 2*Xb*Zd + Zd^2
  return(Xb2)
}

#' Compute E[y - N/2]
compute_kappa <- function(fit, kidx=NULL){
  kappa <- .get_y(fit, kidx) - 0.5 * fit$data$N
  return(kappa)
}

#' Compute E[w] where w are the PG random variables in the augmented model
compute_omega <- function(fit, kidx=NULL){
  omega <- pg_mean(.get_N(fit, kidx), .get_xi(fit, kidx))
  return(omega)
}

#' Compute a scaled posterior mean
#' TODO: add support for non-zero prior mean
#' @param fit a SER object
compute_nu.binser <- function(fit, idx=NULL, kidx=NULL, shift=0){
  mu0 <- .get_mu0(fit, idx=idx)
  tau0 <- 1./.get_var0(fit, idx=idx)

  kappa <- compute_kappa(fit, kidx)
  Zd <- compute_Zd.binser(fit, idx, shift)
  omega <- compute_omega(fit)
  tmp <- kappa - omega * Zd
  nu <- drop(tmp %*% fit$data$X) + (mu0 * tau0)
  return(nu)
}

#' Compute partial posterior variance
#' use this immediately after updating xi
compute_tau <- function(fit){
  omega <- compute_omega(fit)
  tau <- drop(omega %*% fit$data$X^2)
  return(tau)
}


#' update for variational parameter parameter xi q(w) = PG(N, xi)
update_xi.binser <- function(fit){
  Xb2 <- compute_Xb2.binser(fit)
  xi <- sqrt(abs(Xb2))
  return(xi)
}

#' update intercept and fixed effect covariates
update_delta.binser <- function(fit, idx=NULL, kidx=NULL, shift=0){
  Z <- fit$data$Z
  omega <- compute_omega(fit)
  Xb <- fit$data$X %*% Ebeta.binser(fit, idx) + shift
  kappa <- compute_kappa(fit, kidx)
  delta <- solve((omega * t(Z)) %*% Z, t(Z) %*% (kappa - omega*Xb))
  return(delta)
}

#' update SER variational parameters
#' @param nu
#' @param tau
#' @param pi prior probability vector (length p)
.update_b_ser <- function(nu, tau, pi){
  logits <- log(pi) - 0.5 * log(tau) + 0.5 * nu^2/tau
  logits <- logits - logSumExp(logits)
  alpha <- exp(logits)

  post <- list(
    mu = nu/tau,
    var = 1/tau,
    alpha = alpha
  )
  return(post)
}

#' update q(b)
update_b.binser <- function(fit, idx=NULL, kidx=NULL, shift=0){
  nu <- compute_nu.binser(fit, idx, kidx, shift)
  tau0 <- 1 / .get_var0(fit, idx)
  tau <- tau0 + fit$params$tau
  post <- .update_b_ser(nu, tau, .get_pi(fit, idx))
  return(post)
}


#' Variational M step for updating sigma0
#' Optimize ELBO wrt to prior effect standard deviation
update_sigma0.binser <- function(fit, idx=NULL){
  b2 <- .get_alpha(fit, idx) * (.get_mu(fit, idx)^2 + .get_var(fit, idx))
  b <- Ebeta.binser(fit, idx)
  mu0 <- .get_mu0(fit, idx)
  sigma0 <- sqrt(sum(b2 - 2*b*mu0) + mu0^2)
  return(sigma0)
}


#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
#' @returns a vector of length N with the lower bound for each data point
jj_bound.binser <- function(fit, kidx=NULL){
  xi <- .get_xi(fit,)
  Xb <- compute_Xb.binser(fit)
  kappa <- compute_kappa(fit, kidx)
  n <- fit$data$N

  Xb2 <- compute_Xb2.binser(fit)
  omega <- compute_omega(fit)

  bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi) #+ 0.5 * omega * (xi^2 - Xb2)
  return(bound)
}

explicit_elbo.binser <- function(fit){
  Xb <- compute_Xb.binser(fit)
  kappa <- compute_kappa(fit)
  xi <- .get_xi(fit)
  n <- fit$data$N

  Xb2 <- compute_Xb2.binser(fit)
  omega <- compute_omega(fit)

  ll <- -n * log(2) + kappa * Xb - 0.5*Xb2*omega
  kl <- pg_kl(n, xi)

  return(ll - kl)
}

#' More explicit, only works for Bernoulli observations
jj_bound.logistic <- function(fit, kidx=NULL){
  xi <- .get_xi(fit, kidx)
  Xb <- compute_Xb.binser(fit)
  kappa <- compute_kappa(fit, kidx)

  Xb2 <- compute_Xb2.binser(fit)
  omega <- compute_omega(fit)
  bound <- kappa * Xb - 0.5 * xi - 0.25 / xi * tanh(xi/2) * (Xb2 - xi^2) - log(1 + exp(-xi))
  return(bound)
}


compute_elbo.binser <- function(fit){
  jj <- sum(jj_bound.binser(fit))  # E[p(y | w, b)] - KL[q(b) || p(b)]
  kl <- compute_kl.binser(fit)  # KL[q(b) || p(b)]
  return(jj-kl)
}


#' initialize SER
#' @param data a list containing X, Z, y
#' @param mu0 prior effect mean
#' @param simga0 prior effect standard deviation
init.binser <- function(data){
  n <- nrow(data$X)
  p <- ncol(data$X)  # number of covariates to select
  p2 <- ncol(data$Z)  # number of fixed effects + intercept
  # initialize variational parameters
  params <- list(
    alpha = rep(1, p)/p,
    mu = rep(0, p),
    var = rep(1, p),
    delta = rep(0, p2),
    xi = rep(1e-3, n),
    tau = 0 # initialized at end of this function
  )

  # initialize hyper-parameters
  hypers <- list(
    pi = rep(1, p)/p,  # categorical probability of nonzero effect
    mu0 = 0,  # prior mean of effect
    sigma0 = 1,  # prior variance of effect, TODO: include as init arg
    shift = 0
  )

  # TODO: check that data has y, N, X, Z
  fit.init <- list(
    data = data,
    params = params,
    hypers = hypers,
    elbos = c(-Inf)
  )
  fit.init$params$tau <- compute_tau(fit.init)
  return(fit.init)
}

iter.binser <- function(fit, fit_intercept=TRUE, fit_prior_variance=TRUE){
  # update b
  post <- update_b.binser(fit)
  fit$params$mu <- post$mu
  fit$params$var <- post$var
  fit$params$alpha <- post$alpha

  # update intercept/fixed effect covariates
  if(fit_intercept){
    fit$params$delta <- update_delta.binser(fit)
  }
  if(fit_prior_variance){
    fit$hypers$sigma0 <- update_sigma0.binser(fit)
  }

  # update xi
  fit$params$xi <- update_xi.binser(fit)
  fit$params$tau <- compute_tau(fit)

  return(fit)
}

#' Fit the binomial single effect regression
fit.binser <- function(
    data=NULL, fit=NULL,
    maxiter=10,
    tol=1e-3,
    fit_intercept=TRUE,
    fit_prior_variance=TRUE,
    shift = NULL){

  if(is.null(shift)){
    shift <- fit$hypers$shift
  }

  if(is.null(fit)){
    fit <- init.binser(data)
  }

  for(i in 1:maxiter){
    fit <- iter.binser(fit)

    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
    if (.converged(fit, tol)){
      break
    }
  }
  return(fit)
}


