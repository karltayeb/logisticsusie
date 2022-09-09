

#' @title Initialize cEBMF object
#' @description  Initialize cEBMF object
#' @param Y numerical matrix size NxP
#' @param X_l matrix of size NxJ containning covariates affecting the factors
#' @param X_f matrix of size PxT containning covariates affecting the factors
#' @param K numeric number of factors
#'  @param type_noise specify which kind of noise structure is expected, currently three choices. Whether noise constant accross column ('column_wise'), constant 'constant' or constant across rown 'row_wise'
#' @return a cEBMF object
init_cEBMF <- function(Y, X_l,X_f,K=1, type_noise='column_wise' )
{
  cEBMF.obj <- list(
    Y          = Y,
    Y_fit      = 0*Y,
    loading    = matrix(1,nrow=nrow(Y), ncol=K), #first moment
    factor     = matrix(1,nrow=ncol(Y), ncol=K),
    loading2   = matrix(1,nrow=nrow(Y), ncol=K),#second moment
    factor2    = matrix( 1,nrow=ncol(Y), ncol=K),
    tau        = matrix(1, ncol = ncol(Y), nrow=nrow(Y)),
    X_l        = X_l,
    X_f        = X_f,
    K          = K,
    type_noise = type_noise
  )


  class(cEBMF.obj) <- "cEBMF"
  #initialize tau using observed variance
  tau <- update_tau.cEBMF(cEBMF.obj)$tau
  return(cEBMF.obj)
}

#'@param cEBMF a cEBMF object
#'@param k the factor to be updated
#'
#'
update_cEBMF <-  function(cEBMF.obj, k)
{

  #loading update
  l_k <- cal_expected_loading( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                              se      = l_k$s_i,
                              X       = cEBMF.obj$X_l )
  t_fit <- fit.mococomo(t_data)
  fitted_loading <- post_mean_sd.mococomo (t_fit )
  cEBMF.obj$loading[,k] <-  fitted_loading$mean
  cEBMF.obj$loading2[,k] <- fitted_loading$sd^2+ fitted_loading$mean^2


  #factor update

  f_k <- cal_expected_factor( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = f_k$f_j_hat,
                              se      = f_k$s_j,
                              X       = cEBMF.obj$X_f )
  t_fit <- fit.mococomo(t_data)
  fitted_factor <- post_mean_sd.mococomo (t_fit )
  cEBMF.obj$factor[,k] <-  fitted_factor $mean
  cEBMF.obj$factor2[,k] <-  fitted_factor$sd^2+ fitted_factor$mean^2
  return(cEBMF.obj)
}

#'@param cEBMF a cEBMF object
#'@param k component of interest
#'@return list of two l_i_hat the estimate loading, s_i the estimated standard errors
cal_expected_loading <- function( cEBMF.obj,k){

  f       <- cEBMF.obj$factor[,k]
  mat_up  <- cEBMF.obj$Y*cEBMF.obj$tau *f
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

  l       <- cEBMF.obj$loading[,k]
  mat_up  <- cEBMF.obj$Y*cEBMF.obj$tau *l
  mat_low <- cEBMF.obj$tau *(l^2)
  deno    <- apply(mat_low,2,sum)


  f_j_hat <- apply(mat_up,2,sum)/deno
  s_j     <-  1/sqrt(deno)


  out <- list( f_j_hat = f_j_hat,
               s_j     = s_j)
  return( out)
}

cal_fitted_value.cEBMF <- function(cEBMF.obj)
{
   cEBMF.obj$Y_fit <-   Reduce("+",
                               lapply( 1:cEBMF.obj$K,
                                       function(k)
                                         cEBMF.obj$loading[ ,k]%*%t(cEBMF.obj$factor[ ,k])
                                       )
                                )

  return(cEBMF.obj)
}


cal_expected_residuals.cEBMF <- function(cEBMF.obj)
{
  prod_square_firstmom <-   Reduce("+",
                              lapply( 1:cEBMF.obj$K,
                                      function(k)
                                        (cEBMF.obj$loading[ ,k]^2)%*%t(cEBMF.obj$factor[ ,k]^2)
                                    )
                               )
  prod_sectmom <-   Reduce("+",
                                   lapply( 1:cEBMF.obj$K,
                                           function(k)
                                             (cEBMF.obj$loading2[ ,k] )%*%t(cEBMF.obj$factor2[ ,k] )
                                   )
  )

  R2  <- (cEBMF.obj$Y  -  cEBMF.obj$Y_fit)^2 -  prod_square_firstmom+ prod_sectmom

return(R2)

}

#'@title Estimate noise value
#'@'description Estimate noise value for different type of structure
#'@param cEBMF an CEBMF object


update_tau.cEBMF <- function(cEBMF.obj )
{
  R2 <- cal_expected_residuals.cEBMF(cEBMF.obj)
  if(cEBMF.obj$type_noise== 'constant')
  {
    cEBMF.obj$tau <- matrix(1/mean(R2), ncol=ncol(Y), nrow=nrow(Y))
  }
  if(cEBMF.obj$type_noise== 'row_wise' ){
    cEBMF.obj$tau <-  matrix(
                             rep(
                                  1/apply(R2, 1, mean),
                                  each=ncol(Y)
                                  ),
                             ncol= ncol(Y),
                             byrow=TRUE
                             )
  }
  if(cEBMF.obj$type_noise== 'column_wise' ){
    cEBMF.obj$tau <-  matrix(
                               rep(
                                   1/apply(R2, 2, mean),
                                   each=nrow(Y)
                                    ),
                              ncol= ncol(Y),
                              byrow=FALSE
                             )
  }
  return(cEBMF.obj)
}
