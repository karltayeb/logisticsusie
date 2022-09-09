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
#' data <- set_data_mococomo(betahat = sim$betahat,
#'                                se = sim$se ,
#'                                 X = sim$X)
#' #fit mococomo model
#' fit <- fit.mococomo(data, maxiter=20)
#' plot(fit$elbo)
#' .monotone(fit$elbo)
#'
#' #get posterior quantitoes
#' est<- post_mean_sd.mococomo (fit)
#' head(est)
#'  plot( est$mean, data$betahat)
#'
#' #comparison with ash
#'
#' t_ash <- ash(sim $betahat, sim $se, mixcompdist = "normal")
#' post_mean_ash <- t_ash$result$PosteriorMean
#' plot(est$mean, post_mean_ash)
#' # TODO make a more convincing exemple

fit.mococomo <- function(data, maxiter = 100, tol = 1e-3,max_class,  mult = 2) {

  if(!(class(data)=="data_mococomo"))
  {stop("Please provide object of class data_mococomo")}
  fit <- init.mococomo(data,max_class, gridmult )
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


