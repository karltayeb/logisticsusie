#' Get N x K-1 matrix of PG variational parameters
#' TODO: deprecate this when we call update functions directly rather than
#' passing to logistic susie subroutine
get_xi.mnsusie <- function(fit) {
  xi <- do.call(cbind, purrr::map(fit$logreg_list, ~ purrr::pluck(.x, "xi"))) # N x K-1
  return(xi)
}

compute_Xb.mnsusie <- function(fit, data) {
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, ~compute_Xb(.x, data))) # N x K-1
  return(Xb)
}

compute_psi.mnsusie <- function(fit, data) {
  psi <- do.call(cbind, purrr::map(fit$logreg_list, ~compute_psi(.x, data))) # N x K-1
  return(psi)
}

#' Compute Log Prior Assignment Probabilities
#' For each data point return covariate-dependent prior mixture probabilities
#' @param xb a K-1 vector of predictions
#' @return an n x K matrix of log prior probabilities for each data point
.compute_prior_assignment <- function(xb) {
  .predict2logpi(xb)
}

#' Compute Log Prior Assignment Probabilities
#' For each data point return covariate-dependent prior mixture probabilities
#' @param fit a multinomial SuSiE fit object (or MoCoCoMo)
#' @return an n x K matrix of log prior probabilities for each data point
compute_prior_assignment <- function(fit) {
  # TODO alias compute_Xb with predict so that it works with other functions?
  # TODO make sure GLM predict outputs log-odds scale?
  Xb <- compute_Xb.mnsusie(fit) # N x K-1
  res <- do.call(rbind, apply(Xb, 1, .predict2logpi, simplify = F)) # N x K
  return(res)
}


predict_assignment_mnsusie <- function(fit, data) {
  K <- fit$K
  Xb <- compute_psi(fit, data) # N x K-1
  Xi <- get_xi.mnsusie(fit)

  f <- function(xi, xb) {
    tmp <- cumsum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb) + xb
    tmpK <- sum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb)
    jj <- c(tmp, tmpK)
    assignment <- softmax(jj)
    return(assignment)
  }

  assignment <- do.call(rbind, purrr::map(seq(nrow(Xb)), ~ f(Xi[.x, ], Xb[.x, ])))

  return(assignment)
}


compute_elbo.mnsusie <- function(fit, data) {
  elbo <- sum(purrr::map_dbl(fit$logreg_list, ~tail(.x$elbo, 1)))
  return(elbo)
}

update_model.mnsusie <- function(fit, data,
                                 fit_intercept=T,
                                 fit_prior_variance=T,
                                 track_elbo=T){
  for (k in seq(fit$K - 1)) {
    logreg <- fit$logreg_list[[k]]

    # set data for this regression problem
    data$y <- data$Y[, k]
    data$N <- data$Nk[, k]

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


initialize_sbmn_susie <- function(K, n, p, p2, L, mu0=0, var0=1){
  logreg_list <- purrr::map(1:(K-1), ~initialize_binsusie(L, n, p, p2, mu0, var0))

  fit <- list(
    logreg_list = logreg_list,
    K = K, L = L,
    elbo = c(-Inf)
  )
  class(fit) <- 'mnsusie'
  return(fit)
}

data_initialize_sbmn_susie <- function(data, L = 5, mu0=0, var0=1) {
  K <- ncol(data$Y)
  p <- ncol(data$X)
  n <- nrow(data$X)
  p2 <- ncol(data$Z)
  fit <- initialize_sbmn_susie(K, n, p, p2, L, mu0, var0)
}


