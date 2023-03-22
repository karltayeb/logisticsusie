# Categorical from binary

#' Get N x K-1 matrix of PG variational parameters
#' TODO: deprecate this when we call update functions directly rather than
#' passing to logistic susie subroutine
get_xi.cbsusie <- function(fit) {
  xi <- do.call(cbind, purrr::map(fit$logreg_list, ~ purrr::pluck(.x, "xi"))) # N x K-1
  return(xi)
}

compute_Xb.cbsusie <- function(fit, data) {
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, ~compute_Xb(.x, data))) # N x K-1
  return(Xb)
}

compute_psi.cbsusie <- function(fit, data) {
  psi <- do.call(cbind, purrr::map(fit$logreg_list, ~compute_psi(.x, data))) # N x K-1
  return(psi)
}

#' Compute E[log p(z | X, B)]
#' @export
predict_log_assignment_cbsusie <- function(fit, data) {
  K <- fit$K
  Xb <- compute_psi(fit, data) # N x K-1
  Xi <- get_xi.cbsusie(fit)

  f <- function(xi, psi) {
    jj1 <- log(sigmoid(xi)) - 0.5 * xi + 0.5 * psi
    jj0 <- log(sigmoid(xi)) - 0.5 * xi - 0.5 * psi

    # lower bound log[\sigma(\psi_k) \prod_{j \neq k} \sigma(- \psi_j)]
    tmp <- sum(jj0) - jj0 + jj1
    assignment <- tmp - logSumExp(tmp)
    return(assignment)
  }

  assignment <- do.call(rbind, purrr::map(seq(nrow(Xb)), ~ f(Xi[.x, ], Xb[.x, ])))
  return(assignment)
}

#' @export
compute_elbo.cbsusie <- function(fit, data) {
  elbo <- sum(purrr::map_dbl(fit$logreg_list, ~tail(.x$elbo, 1)))
  return(elbo)
}

#' @export
update_model.cbsusie <- function(fit, data,
                                 fit_intercept=T,
                                 fit_prior_variance=T,
                                 track_elbo=T){
  for (k in seq(fit$K - 1)) {
    logreg <- fit$logreg_list[[k]]

    # set data for this regression problem
    data$y <- data$Y[, k]
    data$N <- rep(1, length(data$y))

    # update logreg
    # note we calculate the elbo here since the multinomial susie ELBO is just
    # the sum of it's stick-breaking components
    fit$logreg_list[[k]] <- update_model(fit$logreg_list[[k]], data,
                                         fit_intercept = fit_intercept,
                                         fit_prior_variance = fit_prior_variance,
                                         track_elbo = T
    )
  }

  if(track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}


#' @export
initialize_cbsusie <- function(K, n, p, p2, L, mu0=0, var0=1){
  logreg_list <- purrr::map(1:K, ~initialize_binsusie(L, n, p, p2, mu0, var0))

  fit <- list(
    logreg_list = logreg_list,
    K = K, L = L,
    elbo = c(-Inf)
  )
  class(fit) <- 'cbsusie'
  return(fit)
}

#' @export
data_initialize_cbsusie <- function(data, L = 5, mu0=0, var0=1) {
  K <- ncol(data$Y)
  p <- ncol(data$X)
  n <- nrow(data$X)
  p2 <- ncol(data$Z)
  fit <- initialize_cbsusie(K, n, p, p2, L, mu0, var0)
}


