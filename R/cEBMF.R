
#' @title covariate moderated Empirical Bayes Matrix Factorization
#' @description covariate moderated Empirical Bayes Matrix Factorization
#' @param Y numerical matrix size NxP
#' @param X_l matrix of size NxJ containning covariates affecting the factors
#' @param X_f matrix of size PxT containning covariates affecting the factors
#' @param K numeric number of factors
#' @param type_noise specify which kind of noise structure is expected, currently three choices. Whether noise constant accross column ('column_wise'), constant 'constant' or constant across rown 'row_wise'
#' @param maxit maximum nuber of iteratiob
#' @return a cEBMF object
#'
#'
#'

cEBMF <- function( Y, X_l,X_f,K=1, type_noise='constant',init_type="udv_si", maxit=100  ){

  cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

  for (i in 1:maxit) {
    cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
    if (.converged(fit, tol)) {
      break
    }
  }
  cEBMF.obj <- out_prep.cEBMF(cEBMF.obj)
  return(cEBMF.obj)
}
