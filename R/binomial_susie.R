####
# Binomial SuSiE
###

# Getters------

#' @export
get_alpha.binsusie <- function(fit){
  do.call(rbind, purrr::map(1:fit$L, ~get_alpha(fit$sers[[.x]])))
}

#' @export
get_mu.binsusie <- function(fit){
  do.call(rbind, purrr::map(1:fit$L, ~get_mu(fit$sers[[.x]])))
}

#' @export
get_var.binsusie <- function(fit){
  do.call(rbind, purrr::map(1:fit$L, ~get_var(fit$sers[[.x]])))
}

#' @export
get_delta.binsusie <- function(fit){
  do.call(rbind, purrr::map(1:fit$L, ~get_delta(fit$sers[[.x]])))
}

#' @export
get_var0.binsusie <- function(fit){
  do.call(rbind, purrr::map(1:fit$L, ~get_var0(fit$sers[[.x]])))
}

#' @export
get_mu0.binsusie <- function(fit){
  do.call(rbind, purrr::map(1:fit$L, ~get_mu0(fit$sers[[.x]])))
}

# Computations ----

#' Get Binomial SER coefficients
#' Extract regression coefficients from binomial SER fit
#' @param fit Binomial SER object
#' @return Return E[\beta]
coef.binsusie <- function(fit, data) {
  b <- colSums(get_alpha(fit) * get_mu(fit))
  return(b)
}

#' Expected linear prediction
compute_Xb.binsusie <- function(fit, data) {
  Xb <- Matrix::drop(data$X %*% colSums(get_alpha(fit) * get_mu(fit)))
  return(Xb)
}

#' Expected linear prediction
#' @export
compute_psi.binsusie <- function(fit, data) {
  Xb <- compute_Xb(fit, data)
  Zd <- Matrix::drop(data$Z %*% colSums(get_delta(fit)))
  return(Xb + Zd)
}

#' Second moment of linear prediction
#' @export
compute_Xb2.binsusie <- function(fit, data) {
  B <- get_alpha(fit) * get_mu(fit)
  XB <- data$X %*% t(B)
  Xb <- rowSums(XB)

  B2 <- get_alpha(fit) * (get_mu(fit)^2 + get_var(fit))
  b2 <- colSums(B2)

  Xb2 <- (data$X2 %*% b2)[, 1] + Xb**2 - rowSums(XB^2)
  return(Xb2)
}

#' @export
compute_psi2.binsusie <- function(fit, data){
  Xb2 <- compute_Xb2(fit, data)
  Xb <- compute_Xb(fit, data)
  Zd <- drop(data$Z %*% colSums(get_delta(fit)))
  psi2 <- drop(Xb2 + (2 * Xb * Zd) + Zd^2)
  return(psi2)
}

#' Compute KL(q(\beta) || p(\beta)) for Sum of SERs
#' Note this does not include KL(q(\omega) || p(\omega))
#' since that gets computes as part of the JJ bound
compute_kl.binsusie <- function(fit) {
  kl <- sum(purrr::map_dbl(seq(fit$L), ~ compute_kl(fit$sers[[.x]])))
  return(kl)
}

#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
compute_jj.binsusie <- function(fit, data) {
  xi <- fit$xi
  psi <- compute_psi(fit, data)
  kappa <- data$y - 0.5 * data$N
  n <- data$N

  bound <- n * log(sigmoid(xi)) + (kappa * psi) - (0.5 * n * xi)
  return(bound)
}


compute_elbo.binsusie <- function(fit, data) {
  jj <- sum(compute_jj(fit, data)) # E[p(y | w, b)] - KL[q(b) || p(b)]
  kl <- compute_kl(fit) # KL[q(b) || p(b)]
  return(jj - kl)
}

#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
compute_jj2.binsusie <- function(fit, data) {
  jj <- jj_bound.binsusie(fit, data)
  psi2 <- compute_psi2(fit, data)
  omega <- pg_kl(data$N, fit$xi)
  bound <- jj + 0.5 * omega * (xi^2 - psi2)
  return(bound)
}

compute_elbo2.binsusie <- function(fit, data) {
  jj <- sum(jj_bound2.binsusie(fit, data)) # E[log p(y, w | b) - logp(w)]
  kl <- compute_kl.binsusie(fit) # KL(q(b) || p(b))
  return(jj - kl)
}

###
# Updates ----
###

add_re <- function(mu, mu2, psi, psi2){
  new_mu = mu + psi
  new_mu2 = mu2 + psi2 + (2 * mu * psi)
  return(list(mu=new_mu, mu2=new_mu2))
}

sub_re <- function(mu, mu2, psi, psi2){
  new_mu = mu - psi
  new_mu2 = mu2 - psi2 - (2 * (mu-psi) * psi)
  return(list(mu=new_mu, mu2=new_mu2))
}

#' update alpha, mu, and var
#' @export
update_model.binsusie <- function(fit, data,
                               fit_intercept = TRUE,
                               fit_prior_variance = TRUE,
                               fit_alpha = TRUE,
                               fit_xi = TRUE,
                               update_idx = NULL,
                               track_elbo=TRUE) {
  if (is.null(update_idx)) {
    update_idx <- seq(fit$L)
  }

  re = list(
    mu = compute_psi(fit, data),
    mu2 = compute_psi2(fit, data)
  )
  for (l in update_idx) {
    # remove current effect estimate
    re = sub_re(re$mu,
                re$mu2,
                compute_psi(fit$sers[[l]], data),
                compute_psi2(fit$sers[[l]], data))

    data$shift <- re$mu
    data$shift_var <- re$mu2 - re$mu^2

    # update SER
    fit$sers[[l]]$xi <- fit$xi
    fit$sers[[l]]$tau <- fit$tau
    fit$sers[[l]] <- update_model(fit$sers[[l]], data,
                                  fit_intercept = fit_intercept,
                                  fit_prior_variance = fit_prior_variance,
                                  fit_alpha = fit_alpha,
                                  fit_xi = T)


    # add current effect estimate
    re = add_re(re$mu,
                re$mu2,
                compute_psi(fit$sers[[l]], data),
                compute_psi2(fit$sers[[l]], data))

    fit$xi <- fit$sers[[l]]$xi
    fit$tau <- fit$sers[[l]]$tau

    fit$sers[[l]]$xi <- NULL
    fit$sers[[l]]$tau <- NULL
  }

  if(track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}

# Initialization -------
#' @export
initialize_binsusie <- function(L, n, p, p2, mu0=0, var0=1, pi = rep(1/p, p)){
  if(length(mu0) < L){
    mu0 <- rep(mu0, L)
  }
  if(length(var0) < L){
    var0 <- rep(var0, L)
  }
  sers <- purrr::map(1:L, ~initialize_binser(n, p, p2, mu0[.x], var0[.x], pi))
  xi <- rep(1e-3, n)
  tau <- rep(1, p)
  fit <- list(sers=sers, xi=xi, tau=tau, L=L, elbo=-Inf)
  class(fit) <- 'binsusie'
  return(fit)
}

#' @export
data_initialize_binsusie <- function(data, L, mu0=0, var0=1){
  n <- nrow(data$X)
  p <- ncol(data$X)
  p2 <- ncol(data$Z)
  fit <- initialize_binsusie(L, n, p, p2, mu0, var0)
  fit$xi <- update_xi(fit, data)
  fit$tau <- compute_tau(fit, data)
  return(fit)
}

#' BinSuSiE Wrapup
#' Compute fit summaries (PIPs, CSs, etc) and
#' Organize binsusie fit so that it is compatable with `susieR` function
binsusie_wrapup <- function(fit, prior_tol = 0) {
  class(fit) <- c("binsusie", "susie")

  # TODO put back into original X scale
  X <- fit$data$X

  fit$alpha <- fit$params$alpha
  fit$mu <- fit$params$mu
  fit$mu2 <- fit$params$mu^2 + fit$params$var
  fit$V <- fit$hypers$prior_variance

  colnames(fit$alpha) <- colnames(X)
  colnames(fit$mu) <- colnames(X)

  fit$pip <- susieR::susie_get_pip(fit, prior_tol = prior_tol)
  names(fit$pip) <- colnames(X)

  fit$sets <- compute_cs(fit$alpha)
  fit$intercept <- colSums(fit$params$delta)[1]
  return(fit)
}

# help -----

#' @export
binsusie_get_pip <- function(fit, prune_by_cs = FALSE, prior_tol = 1e-09) {
  # TODO: filter out components that don't need to be included (see Gao's suggestion)
  include_idx <- which(fit$hypers$prior_variance > prior_tol)
  if (length(include_idx) > 0) {
    pip <- susieR::susie_get_pip(fit$params$alpha[include_idx, , drop = F])
  } else {
    pip <- rep(0, dim(fit$params$alpha)[2])
  }
  return(pip)
}

#' @export
binsusie_plot <- function(fit, y = "PIP") {
  res <- with(fit, list(alpha = params$alpha, pip = pip, sets = sets))
  class(res) <- "susie"
  susieR::susie_plot(res, y)
}

#' Binomial SuSiE
#' Fit Binomial SuSiE via coordinate ascent variational inference
#' @param X a n x p matrix of covariates
#' @param y an n vector of integer counts, bernoulli/binomial observations
#' @param N the number of binomial trials, defaults to 1, may be a scalar or vector of length n
#' @param Z fixed effect covaraites (including intercept). If null just a n x 1 matrix of ones
#' @param scale if TRUE, scale the columns of X to unit variate
#' @param center if TRUE, center the columns of X to mean zero
#' @param prior_mean the prior mean of each non-zero element of b. Either a scalar or vector of length L.
#' @param prior_variance the prior variance of each non-zero element of b. Either a scalar or vector of length L. If `estimate_prior_variance=TRUE` the value provides an initial estimate of prior variances
#' @param prior_weights prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix
#' @param intercept
#' @param estimate_prior_variance
#' @param s_init a logistic susie object to initialize with, NOTE if non-null, we ignore `prior_mean`, `prior_variance`, and `prior_weights`
#' @param returns a fit Binomial SuSiE model, which is compatable with summary functions from `susieR` package
#' @export
binsusie <- function(X,
                     y,
                     N = rep(1, length(y)), # number of trials for each
                     Z = NULL,
                     L = min(10, ncol(X)),
                     scale = FALSE,
                     center = TRUE,
                     prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
                     prior_variance = 1.0, # TODO: make scaled prior variance?
                     prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect... = hypers$pi
                     # null_weight = NULL,
                     # standardize = TRUE,
                     intercept = TRUE,
                     estimate_prior_variance = TRUE,
                     # estimate_prior_method = c("optim", "EM", "simple"),  right now only EM is implimented
                     check_null_threshold = 0,
                     prior_tol = 1e-09,
                     prune = FALSE,
                     s_init = NULL, # previous fit with which to initialize NO
                     coverage = 0.95,
                     min_abs_corr = 0.5,
                     max_iter = 100,
                     tol = 0.001,
                     verbose = FALSE,
                     # track_fit = FALSE,
                     # refine = FALSE,
                     n_purity = 100) {

  data <- binsusie_prep_data(X, y, N, Z, center=center, scale=scale)
  fit <- data_initialize_binsusie(data, L, mu0 = prior_mean, var0 = prior_variance)
  fit <- fit_model(fit, data, max_iter=max_iter, tol=tol)

  fit$alpha <- do.call(rbind, purrr::map(fit$sers, ~.x$alpha))
  fit$mu <- do.call(rbind, purrr::map(fit$sers, ~.x$mu))
  fit$var <- do.call(rbind, purrr::map(fit$sers, ~.x$var))

  # model pruning-- removes irrelevant components
  # if (prune) {
  #   fit <- prune_model(fit, check_null_threshold, intercept, estimate_prior_variance, tol = tol)
  # }
  # wrapup (computing PIPs, CSs, etc)
  # fit <- binsusie_wrapup(fit, prior_tol)
  return(fit)
}
