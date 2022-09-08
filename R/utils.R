

#' @title Prepraring data for mococomo fit
#' @details Prepraring data for mococomo fit
#'
#' @param betahat numeric, vecotr of regression coefficients  list of containing the following element see
#' @param maxiter numeric, maximum numerous of itration set to 100 by defaults
#' @param tol tolreance in term of chjange in ELBO value for stopping criterion
#' @export
#' @example
#' see \link{\code{fit.mococomo}}
set_data_mococomo <- function (beta, se, X, ... ){

  if( !is.numeric( beta))
  {
    stop("beta should be numercial vector")
  }
  if( !is.numeric( se))
  {
    stop("se should be numercial vector")
  }
  if( !is.numeric( X))
  {
    stop("X should be numercial vector")
  }
  if( !is.matrix( beta))
  {
    stop("Xshould be a matrix")
  }
  if( ! (sum(c(length(se)==length(betahat), length(se)==nrow(X)))==2  ))
  {
    stop(" The number of lines in X should be equal to the number of entries in Betahat and se ")
    }


  dat <-  list(betahat= beta,
               se=se ,
               X=X)
  class(data) <- "data_mococomo"
}
