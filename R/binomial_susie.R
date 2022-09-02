####
# Binomial SuSiE
###


###
# Computations
##

#' Expected linear prediction
compute_Xb.binsusie <- function(fit){
  Xb <- drop(fit$data$X %*% colSums(fit$params$alpha * fit$params$mu))
  Zd <- drop(fit$data$Z %*% colSums(fit$params$delta))
  return(Xb + Zd)
}


#' Second moment of linear prediction
compute_Xb2.binsusie <- function(fit){
  B <- .get_alpha(fit) * .get_mu(fit)
  XB <- fit$data$X %*% t(B)
  Xb <- rowSums(XB)
  Zd <- drop(fit$data$Z %*% colSums(fit$params$delta))

  B2 <- fit$params$alpha * (fit$params$mu^2 + fit$params$var)
  b2 <- colSums(B2)

  Xb2 <- fit$data$X^2 %*% b2 + Xb**2 - rowSums(XB^2)
  Xb2 <- drop(Xb2 + 2*Xb*Zd + Zd^2)
  return(Xb2)
}

#' Compute KL(q(\beta) || p(\beta)) for Sum of SERs
#' Note this does not include KL(q(\omega) || p(\omega))
#' since that gets computes as part of the JJ bound
compute_kl.binsusie <- function(fit){
  kl <- sum(purrr::map_dbl(seq(fit$hypers$L), ~ .compute_ser_kl(
    .get_alpha(fit, .x),
    .get_pi(fit, .x),
    .get_mu(fit, .x),
    .get_var(fit, .x),
    .get_mu0(fit, .x),
    .get_var0(fit, .x)
  )))
  return(kl)
}


#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
jj_bound.binsusie <- function(fit){
  xi <- fit$params$xi
  Xb <- compute_Xb.binsusie(fit)
  kappa <- compute_kappa(fit)
  n <- fit$data$N

  Xb2 <- compute_Xb2.binsusie(fit)
  omega <- compute_omega(fit)

  bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi) #+ 0.5 * omega * (xi^2 - Xb2)
  return(bound)
}


compute_elbo.binsusie <- function(fit){
  jj <- sum(jj_bound.binsusie(fit))
  kl <- compute_kl.binsusie(fit)
  return(jj-kl)
}

###
# Updates
###


#' update alpha, mu, and var
update_b.binsusie <- function(fit, fit_intercept=TRUE, fit_prior_variance=TRUE){
  shift <- compute_Xb.binsusie(fit)
  for(l in seq(fit$hypers$L)){
    # remove current effect estimate
    shift <- shift - compute_Xb.binser(fit, idx=l)

    # update SER
    post_l <- update_b.binser(fit, idx=l, shift=shift)
    fit$params$mu[l,] <- post_l$mu
    fit$params$var[l,] <- post_l$var
    fit$params$alpha[l,] <- post_l$alpha

    # update intercept/fixed effect covariates
    if(fit_intercept){
      fit$params$delta[l,] <- update_delta.binser(fit, idx=l, shift=shift)
    }

    # update sigma0
    if(fit_prior_variance){
      fit$hypers$sigma0[l] <- update_sigma0.binser(fit, idx=l)
    }

    # add current effect estimate
    shift <- shift + compute_Xb.binser(fit, idx=l)
  }
  return(fit)
}


#' update for variational parameter parameter xi q(w) = PG(N, xi)
update_xi.binsusie <- function(fit){
  Xb2 <- compute_Xb2.binsusie(fit)
  xi <- sqrt(abs(Xb2))
  return(xi)
}

###
# Fitting
###

#' initialize SER
init.binsusie <- function(data, L=5){
  n <- nrow(data$X)
  p <- ncol(data$X)
  p2 <- ncol(data$Z)

  # initialize variational parameters
  params <- list(
    alpha = matrix(rep(1, p*L)/p, nrow=L),
    mu = matrix(rep(0, L*p), nrow=L),
    var = matrix(rep(1, L*p), nrow=L),
    delta = matrix(rep(0, p2*L, nrow=L)),
    xi = rep(1e-3, n),
    tau = 0 # initialized at end of this function
  )

  # initialize hyper-parameters
  hypers <- list(
    L = L,
    pi = matrix(rep(1, p*L)/p, nrow=L),  # categorical probability of nonzero effect
    mu0 = rep(0, L),  # prior mean of effect
    sigma0 = rep(1, L),  # prior variance of effect, TODO: include as init arg,
    shift = rep(0, n)
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


#' Fit the binomial single effect regression
fit.binsusie <- function(
    data,
    maxiter=10,
    tol=1e-3,
    fit_intercept=TRUE,
    fit_prior_variance=TRUE){

  fit <- init.binsusie(data)

  for(i in 1:maxiter){
    # update b
    fit <- update_b.binsusie(fit, fit_intercept, fit_prior_variance)

    # update xi
    fit$params$xi <- update_xi.binsusie(fit)
    fit$params$tau <- compute_tau(fit)

    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.binsusie(fit))
    if (.converged(fit, tol)){
      break
    }
  }
  return(fit)
}



