####
# Binomial SuSiE
###


###
# Computations
##

compute_Xb.binsusie <- function(fit){
  Xb <- drop(fit$data$X %*% colSums(fit$params$alpha * fit$params$mu))
  return(Xb)
}

compute_kl.binsusie <- function(fit){
  kl <- sum(purrr::map_dbl(seq(fit$hypers$L), ~ .compute_binomial_ser_kl(
    .get_alpha(fit, .x),
    .get_pi(fit, .x),
    .get_mu(fit, .x),
    .get_var(fit, .x),
    .get_mu0(fit, .x),
    .get_var0(fit, .x)
  )))
  return(kl)
}

###
# Updates
###


#' update alpha, mu, and var
update_b.binsusie <- function(fit){

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
    delta = matrix(rep(1, p2*L, nrow=L)),
    xi = rep(1e-3, n),
    tau = 0 # initialized at end of this function
  )

  # initialize hyper-parameters
  hypers <- list(
    L = L,
    pi = matrix(rep(1, p*L)/p, nrow=L),  # categorical probability of nonzero effect
    mu = 0,  # prior mean of effect
    sigma0 = 1.,  # prior variance of effect, TODO: include as init arg,
    offset = rep(0, n)
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


#' update one of the single effects
update_ser.binsusie <- function(fit, k){
  # subtract current estimate
  #fit$hypers$offset <<- fit$hypers$offset - drop(fit$data$X %*% (fit$params$alpha[k,] * fit$params$mu[k,]))
  fit$hypers$offset <- fit$hypers$offset - drop(fit$data$X %*% (fit$params$alpha[k,] * fit$params$mu[k,]))

  # update SER, give update_b.binser the offset and current index
  fit$hypers$offset <- offset
  fit$hypers$idx <- k

  post <- update_b.binser(fit)
  fit$params$alpha[k,] <<- post$alpha
  fit$params$var[k,] <<- post$var
  fit$params$mu[k,] <<- post$mu

  # add back current component
  fit$hypers$offset <<- fit$hypers$offset + drop(fit$data$X %*% (fit$params$alpha[k,] * fit$params$mu[k,]))
}

iter.binsusie <- function(fit){

  offset.init <- compute_Xb.binsusie(fit)
  purrr::accumulate(seq(fit$hypers$L), iter.binsusie, .init=offset.init)
}


#' Fit the binomial single effect regression
fit.binsusie <- function(
    data,
    maxiter=10,
    tol=1e-3,
    fit_intercept=TRUE,
    fit_prior_variance=TRUE){

  fit <- init.binser(data)

  for(i in 1:maxiter){
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

    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.binser(fit))
    if (.converged(fit, tol)){
      break
    }
  }
  return(fit)
}



