

#' @title Initialize cEBMF object
#' @description  Initialize cEBMF object
#' @param Y numerical matrix size NxP
#' @param X_l matrix of size NxJ containning covariates affecting the factors
#' @param X_f matrix of size PxT containning covariates affecting the factors
#' @param K numeric number of factors
#' @return a cEBMF object
init_cEBMF <- function(Y, X_l,X_f,K=1)
{
  cEBMF.obj <- list(
    Y=Y,
    loading =  matrix(1,nrow=nrow(Y), ncol=K),
    factor  =  matrix( 1,nrow=ncol(Y), ncol=K),
    tau     = matrix(1, ncol = ncol(Y), nrow=nrow(Y)),
    X_l     = X_l,
    X_f     = X_f
  )
  class(cEBMF.obj) <- "cEBMF"
  return(cEBMF.obj)
}

#'@param cEBMF a cEBMF object
#'@param k the factor to be updated
#'
#'
update_cEBMF <-  function(cEBMF.obj, k)
{
  l_k <- cal_expected_loading( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                              se      = l_k$s_i,
                              X       = cEBMF.obj$X_l )
  t_fit <- fit.mococomo(t_data)
  l_i_fit <- post_mean.mococomo (fit , t_data)

}

#'@param cEBMF a cEBMF object
#'@param k component of interest
#'@return list of two l_i_hat the estimate loading, s_i the estimated standard errors
cal_expected_loading <- function( cEBMF.obj,k){

  f       <- cEBMF.obj$f[,k]
  mat_up  <- cEBMF.obj$Y*cEBMF.obj$tau *cEBMF.obj$f[,k]
  mat_low <- cEBMF.obj$tau *(f^2)
  deno    <- apply(mat_low,1,sum)


  l_i_hat <- apply(mat_up,1,sum)/deno
  s_i     <-  1/sqrt(deno)


  out <- list( l_i_hat = l_i_hat,
               s_i     = s_i)
  return( out)
}



#'@param cEBMF a cEBMF object
#'@param k  component of interest
#'@return list of two l_i_hat the estimate loading, s_i the estimated standard errors
cal_expected_factor <- function( cEBMF.obj,k){

  l       <- cEBMF.obj$l[,k]
  mat_up  <- cEBMF.obj$Y*cEBMF.obj$tau *l
  mat_low <- cEBMF.obj$tau *(l^2)
  deno    <- apply(mat_low,2,sum)


  f_j_hat <- apply(mat_up,2,sum)/deno
  s_j     <-  1/sqrt(deno)


  out <- list( f_j_hat = l_j_hat,
               s_j     = s_j)
  return( out)
}


post_mean.mococomo (fit , t_data)
