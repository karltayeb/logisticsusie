# A common interface for binomial SER, binomial SuSiE and multinomial SuSiE

predict_model <- function(x, ...){
  UseMethod("predict_model", x)
}
predict_model.default <- function(fit, data){
  compute_psi(fit, data)
}

update_model <- function(x, ...){
  UseMethod("update_model", x)
}

#' Fit model
#'
#' An outer loop for `update_model` that monitors convergence, etc.
fit_model <- function(x, ...){
  UseMethod("fit_model", x)
}

fit_model.default <- function(fit, data, ..., max_iter=100, tol=1e-5){
  for(i in 1:max_iter){
    fit <- update_model(fit, data, ...)

    # Check convergence, assumes fit is tracking elbo
    if(abs(diff(tail(fit$elbo, 2))) < tol){
      message('converged')
      break
    }
  }
  return(fit)
}


compute_coef <- function(x, ...){
  UseMethod("compute_coef")
}
compute_coef.default <- function(fit){
  b <- Matrix::drop(get_alpha(fit) * get_mu(fit))
}

compute_Xb <- function(x, ...){
  UseMethod("compute_Xb")
}
compute_psi <- function(x, ...){
  UseMethod("compute_psi")
}

compute_Xb2 <- function(x, ...){
  UseMethod("compute_Xb2")
}
compute_psi2 <- function(x, ...){
  UseMethod("compute_psi2")
}

coef <- function(x, ...){
  UseMethod("coef")
}

get_alpha <- function(x, ...){
  UseMethod("get_alpha", x)
}

get_mu <- function(x, ...){
  UseMethod("get_mu", x)
}

get_pi <- function(x, ...){
  UseMethod("get_pi", x)
}

get_var <- function(x, ...){
  UseMethod("get_var", x)
}

get_delta <- function(x, ...){
  UseMethod("get_delta", x)
}

get_mu0 <- function(x, ...){
  UseMethod("get_mu0", x)
}

get_var0 <- function(x, ...){
  UseMethod("get_var0", x)
}

compute_jj <- function(x, ...){
  UseMethod("compute_jj", x)
}

compute_kl <- function(x, ...){
  UseMethod("compute_kl", x)
}

compute_elbo <- function(x, ...){
  UseMethod("compute_elbo", x)
}
compute_elbo.default <- function(fit, data) {
  jj <- sum(compute_jj(fit, data)) # E[p(y | w, b)] - KL[q(b) || p(b)]
  kl <- compute_kl(fit) # KL[q(b) || p(b)]
  return(jj - kl)
}

### JJ functions

#' Compute E[y - N/2]
compute_kappa <- function(data) {
  kappa <- data$y - 0.5 * data$N
  return(kappa)
}


#' update for variational parameter parameter xi q(w) = PG(N, xi)
update_xi <- function(fit, data) {
  psi2 <- compute_psi2(fit, data)
  xi <- sqrt(abs(psi2))
  return(xi)
}

#' Compute E[w] where w are the PG random variables in the augmented model
compute_omega <- function(fit, data) {
  omega <- pg_mean(data$N, fit$xi)
  return(omega)
}

#' Compute partial posterior variance
#' use this immediately after updating xi
compute_tau <- function(fit, data) {
  omega <- compute_omega(fit, data)
  tau <- Matrix::drop(omega %*% data$X2)
  return(tau)
}

#' KL(q(b, \omega) || p(\beta, \omega))
compute_ser_kl <- function(alpha, pi, mu, var, mu0, var0) {
  kl <- categorical_kl(alpha, pi)
  kl <- kl + sum(alpha * normal_kl(
    mu, var, mu0, var0
  ))
  return(kl)
}
