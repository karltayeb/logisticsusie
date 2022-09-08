# Impelements covariate moderated ASH "MOre COmponents COvariate MOderated"

#' @title Fucntion implementating the mococomo mode
#' @details Fucntion implementating the mococomo mode
#'
#' @param data an object of class data_mococomo  see \link{\code{set_data_mococomo}}
#' @param maxiter numeric, maximum numerous of itration set to 100 by defaults
#' @param tol tolreance in term of chjange in ELBO value for stopping criterion
#' @export
#' @example
#' #Simulate data under the mococomo model
#' sim  <- sim_twococomo()
#' #preparing the data
#' data <-  set_data(betahat = data$betahat,
#'                        se = data$se ,
#'                        X  = data$X)
#' #fit mococomo model
#' fit <- fit.mococomo(data, maxiter=20)
#' plot(fit$elbo)
#' .monotone(fit$elbo)

fit.mococomo <- function(data, maxiter = 100, tol = 1e-3) {

  if(!(class(data)=="data_mococomo"))
  fit <- init.mococomo(data)
  fit$elbo <- compute_elbo.mococomo(fit)
  for (i in 1:maxiter) {
    fit <- iter.mococomo(fit, is.even(i), is.odd(i))
    fit$elbo <- c(fit$elbo, compute_elbo3.mococomo(fit))

    # print(paste('asgn:', is.even(i), 'logreg:', is.odd(i), 'elbo: ', tail(fit$elbo, 1)))
    if (.converged(fit, tol)) {
      break
    }
  }
  return(fit)
}


