#' Get N x K-1 matrix of PG variational parameters
#' TODO: deprecate this when we call update functions directly rather than
#' passing to logistic susie subroutine
get_xi.mnsusie <- function(fit) {
  xi <- do.call(cbind, purrr::map(fit$logreg_list, ~ purrr::pluck(.x, "params", "xi"))) # N x K-1
  return(xi)
}

#' Conditional number of binomial trials
#' @param Y an N x K matrix of assignment probabilities
compute_N.mnsusie <- function(Y) {
  cumY <- do.call(rbind, apply(Y, 1, cumsum, simplify = F))
  K <- dim(Y)[2]
  N <- cumY[, K] - cumY + Y
  return(N)
}

compute_Xb.mnsusie <- function(fit) {
  Xb <- do.call(cbind, purrr::map(fit$logreg_list, compute_Xb.binsusie)) # N x K-1
  return(Xb)
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



compute_jj_bound.mnsusie <- function(fit) {
  K <- length(fit$f_list)
  Xb <- compute_Xb.mnsusie(fit) # N x K-1
  Xi <- get_xi.mnsusie(fit)

  f <- function(xi, xb) {
    tmp <- cumsum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb) + xb
    tmpK <- sum(log(sigmoid(xi)) - 0.5 * xi - 0.5 * xb)
    jj <- c(tmp, tmpK)
    return(jj)
  }
  jj <- do.call(rbind, purrr::map(seq(nrow(Xb)), ~ f(Xi[.x, ], Xb[.x, ])))
  return(jj)
}


compute_elbo.mnsusie <- function(fit, from_init_moco=FALSE) {


  if(from_init_moco){
    K <- length(fit$logreg_list) + 1

    # update omega so jj bound is tight
    for (k in seq(K - 1)) {
      fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit ,k=k)
      fit$logreg_list[[k]]$params$tau <- compute_tau(fit ,
                                                     k=k)
    }

    elbo <- Reduce ("+", lapply( 1:length(fit$logreg_list) ,
                                 function(k) compute_elbo.binsusie(fit,k=k)
    )
    )


  }else{
    K <- length(fit$logreg_list) + 1

    # update omega so jj bound is tight
    for (k in seq(K - 1)) {
      fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit$logreg_list[[k]])
      fit$logreg_list[[k]]$params$tau <- compute_tau(fit$logreg_list[[k]])
    }

    elbo <- sum(purrr::map_dbl(fit$logreg_list, compute_elbo.binsusie))
  }

  return(elbo)
}




init.mnsusie <- function(data, L = 5, from_init_moco=FALSE) {
  data$X2 <- data$X^2

  K <- dim(data$y)[2]
  p <- ncol(data$X)
  n <- nrow(data$X)

  # initialize posterior assignment
  Y <- data$y
  N <- compute_N.mnsusie(Y)

  data$y <- Y
  data$N <- N

  # initialize K-1 logistic SuSiE
  .init_logreg <- function(y, N ) {
    dat <- list(
      X = data$X,
      Z = data$Z,
      y = y,
      N = N
    )
    fit <- init.binsusie(dat,
                         L = L,
                         from_init_moco=from_init_moco)
    return(fit)
  }
  logreg_list <- purrr::map(seq(K - 1), ~ .init_logreg(Y[, .x], N[, .x]))

  if( from_init_moco){
    fit <- list(
      logreg_list = logreg_list,
      elbo = c(-Inf)
    )
  }else{
    fit <- list(
      data = data,
      logreg_list = logreg_list,
      elbo = c(-Inf)
    )
  }


  return(fit)
}


iter.mnsusie <- function(fit,
                         fit_intercept = T,
                         fit_prior_variance = T,
                         from_init_moco = FALSE) {


  if(from_init_moco){
    K <- length(fit$logreg_list) + 1

    logreg_list <- list()
    for (k in seq(K - 1)) {
      # logreg <- fit$logreg_list[[k]]

      # update logreg
      fit$logreg_list[[k]] <- iter.binsusie(fit=fit,
                                            k=k,
                                            fit_intercept = fit_intercept,
                                            fit_prior_variance = fit_prior_variance
      )$logreg_list[[k]]

    }


  }else{
    K <- length(fit$logreg_list) + 1

    logreg_list <- list()
    for (k in seq(K - 1)) {
      logreg <- fit$logreg_list[[k]]

      # update logreg
      logreg_list[[k]] <- iter.binsusie(logreg,
                                        fit_intercept = fit_intercept,
                                        fit_prior_variance = fit_prior_variance
      )
    }


  }


  return(fit)
}


#' Fit the binomial single effect regression
fit.mnsusie <- function(data,
                        L = 5,
                        maxiter = 10,
                        tol = 1e-3,
                        fit_intercept = TRUE,
                        fit_prior_variance = TRUE) {
  fit <- init.mnsusie(data, L = L)

  for (i in 1:maxiter) {
    fit <- iter.mnsusie(fit, fit_intercept, fit_prior_variance)
    # update elbo
    fit$elbo <- c(fit$elbo, compute_elbo.mnsusie(fit))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
}
