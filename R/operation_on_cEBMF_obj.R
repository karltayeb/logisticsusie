#### TODO
#elbo
#Remove factor which prior is null
#greedy and backfitting implementation

#' @title Initialize cEBMF object
#' @description  Initialize cEBMF object
#' @param Y numerical matrix size NxP
#' @param X_l matrix of size NxJ containning covariates affecting the factors
#' @param X_f matrix of size PxT containning covariates affecting the factors
#' @param K numeric number of factors
#'  @param type_noise specify which kind of noise structure is expected, currently three choices. Whether noise constant accross column ('column_wise'), constant 'constant' or constant across rown 'row_wise'
#'  @param init_type specify initialisation method for loading and factor, methods includes: "udv", "udv_si_svd", "si" (softimpute), "ones" and "rand_udv"
#' @return a cEBMF object
init_cEBMF <- function(Y, X_l,X_f,K=1, type_noise='constant',init_type="udv_si" )
{


  if( missing(X_l))
  {
    X_l <- matrix(rnorm(nrow(Y)*10), nrow= nrow(Y))
  }
  if( missing(X_f))
  {
    X_f <- matrix(rnorm(ncol(Y)*10), nrow= ncol(Y))
  }

  if (init_type=="ones"){
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
      type_noise = type_noise,
      elbo       = c(),
      param_loading = list(), #susie param for each loading
      KL_l          = list(),
      param_factor  = list(),  #susie param for each factor
      KL_f          = list()

    )
  }
  if(init_type=="rand_udv"){
    s <- svd( matrix(rnorm(prod(dim(Y))),
                     ncol= ncol(Y)
                    ),
              K,K
             )


    if(K==1){
      l <- (s$u[,1])*(s$d[1])
    }else{
      l <-s$u[,1:K]%*%diag(s$d[1:K])
    }
    f <- s$v[,1:K]


    cEBMF.obj <- list(
      Y          = Y,
      Y_fit      = 0*Y,
      loading    = l  , #first moment
      factor     = f,
      loading2   = l^2,#second moment
      factor2    = f^2,
      tau        = matrix(1, ncol = ncol(Y), nrow=nrow(Y)),
      X_l        = X_l,
      X_f        = X_f,
      K          = K,
      type_noise = type_noise,
      elbo       = c(),
      param_loading = list(), #susie param for each loading
      KL_l          = list(),
      param_factor  = list(),  #susie param for each factor
      KL_f          = list()

    )
  }
  if(init_type=="udv"){

    s <- svd(Y, K,K)
    if(K==1){
      l <- (s$u[,1])*(s$d[1])
    }else{
      l <-s$u[,1:K]%*%diag(s$d[1:K])
    }
    f <- s$v[,1:K]


    cEBMF.obj <- list(
      Y          = Y,
      Y_fit      = 0*Y,
      loading    = l  , #first moment
      factor     = f,
      loading2   = l ^2,#second moment
      factor2    = f ^2,
      tau        = matrix(1, ncol = ncol(Y), nrow=nrow(Y)),
      X_l        = X_l,
      X_f        = X_f,
      K          = K,
      type_noise = type_noise,
      elbo       = c(),
      param_loading = list(), #susie param for each loading
      KL_l          = list(),
      param_factor  = list(),  #susie param for each factor
      KL_f          = list()

    )
  }
  if(init_type=="udv_si"){
    suppressWarnings(
      s <- softImpute(Y, rank.max=K,type="als", lambda=0)
    )
    if(K==1){
      l <- (s$u )*(s$d[1])
      f <- s$v
    }else{
      l <-s$u[,1:K]%*%diag(s$d[1:K])
      f <- s$v[,1:K]
    }


    cEBMF.obj <- list(
      Y          = Y,
      Y_fit      = 0*Y,
      loading    = l , #first moment
      factor     = f,
      loading2   = l^2,#second moment
      factor2    = f^2,
      tau        = matrix(1, ncol = ncol(Y), nrow=nrow(Y)),
      X_l        = X_l,
      X_f        = X_f,
      K          = K,
      type_noise = type_noise,
      elbo       = c(),
      param_loading = list(), #susie param for each loading
      KL_l          = list(),
      param_factor  = list(),  #susie param for each factor
      KL_f          = list()

    )
  }
  if(init_type=="udv_si_svd"){
    suppressWarnings(
      s <- softImpute(Y, rank.max=K,type="svd", lambda=0)
    )
    if(K==1){
      l <- (s$u )*(s$d[1])
      f <- s$v
    }else{
      l <-s$u[,1:K]%*%diag(s$d[1:K])
      f <- s$v[,1:K]
    }

    cEBMF.obj <- list(
      Y          = Y,
      Y_fit      = 0*Y,
      loading    = l , #first moment
      factor     = f,
      loading2   = l^2,#second moment
      factor2    = f^2,
      tau        = matrix(1, ncol = ncol(Y), nrow=nrow(Y)),
      X_l        = X_l,
      X_f        = X_f,
      K          = K,
      type_noise = type_noise,
      elbo       = c(),
      param_loading = list(), #susie param for each loading
      KL_l          = list(),
      param_factor  = list(),  #susie param for each factor
      KL_f          = list()

    )
  }
  class(cEBMF.obj) <- "cEBMF"
  #initialize tau using observed variance
  tau <- update_tau.cEBMF(cEBMF.obj)$tau
  return(cEBMF.obj)
}


#'@param cEBMF a cEBMF object
#'@param k component of interest
#'@return list of two l_i_hat the estimate loading, s_i the estimated standard errors
cal_expected_loading <- function( cEBMF.obj, Rk,k){

  if( cEBMF.obj$K==1){
    f       <- cEBMF.obj$factor
    f2      <- cEBMF.obj$factor2
  }else{
    f       <- cEBMF.obj$factor[,k]
    f2      <- cEBMF.obj$factor2[,k]
  }

  mat_up  <- sweep(Rk*cEBMF.obj$tau ,2,f, "*")
  mat_low <- sweep(cEBMF.obj$tau,2,  f2 , "*")
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
cal_expected_factor <- function( cEBMF.obj, Rk,k){


  if( cEBMF.obj$K==1){
    l       <- cEBMF.obj$loading
    l2      <- cEBMF.obj$loading2
  }else{
    l       <- cEBMF.obj$loading[,k]
    l2      <- cEBMF.obj$loading2[,k]
  }

  mat_up  <- Rk*cEBMF.obj$tau *l
  mat_low <- cEBMF.obj$tau *(l2)
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
  if(cEBMF.obj$K==1){
    prod_square_firstmom <- (cEBMF.obj$loading^2)%*%t(cEBMF.obj$factor^2)

    prod_sectmom <-   (cEBMF.obj$loading2)%*%t(cEBMF.obj$factor2)
  }else{
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
  }



  R2  <- (cEBMF.obj$Y  -  cEBMF.obj$Y_fit)^2 -  prod_square_firstmom+ prod_sectmom

  return(R2)

}


cal_partial_residuals.cEBMF <- function(cEBMF.obj,k)
{
  K <- cEBMF.obj$K
  if (K > 1){

    id_k <- (1:K)[ -k ]
    partial_pred <-   Reduce("+",
                             lapply( id_k,
                                     function(k)
                                       (cEBMF.obj$loading[ ,k])%*%t(cEBMF.obj$factor[ ,k])
                             )

    )
  }else{
    partial_pred <- 0*cEBMF.obj$Y
  }


  Rk  <-  cEBMF.obj$Y  -   partial_pred

  return(Rk)

}


#'@title Estimate noise value
#'@'description Estimate noise value for different type of structure
#'@param cEBMF an CEBMF object


update_tau.cEBMF <- function(cEBMF.obj )
{
  R2 <- cal_expected_residuals.cEBMF(cEBMF.obj)
  if(cEBMF.obj$type_noise== 'constant')
  {
    cEBMF.obj$tau <- matrix(1/mean(R2),
                            ncol=ncol(cEBMF.obj$Y),
                            nrow=nrow(cEBMF.obj$Y)
                            )
  }
  if(cEBMF.obj$type_noise== 'row_wise' ){
    cEBMF.obj$tau <-  matrix(
      rep(
        1/apply(R2, 1, mean),
        each=ncol(cEBMF.obj$Y)
      ),
      ncol= ncol(cEBMF.obj$Y),
      byrow=TRUE
    )
  }
  if(cEBMF.obj$type_noise== 'column_wise' ){
    cEBMF.obj$tau <-  matrix(
      rep(
        1/apply(R2, 2, mean),
        each=nrow(cEBMF.obj$Y)
      ),
      ncol= ncol(cEBMF.obj$Y),
      byrow=FALSE
    )
  }
  return(cEBMF.obj)
}

#'@param cEBMF a cEBMF object

cEBMF_iter <- function(cEBMF.obj){

  for( k in 1:cEBMF.obj$K){

    cEBMF.obj <- update_cEBMF(cEBMF.obj, k)

  }
  cEBMF.obj <- update_tau.cEBMF (cEBMF.obj)
  #cEBMF.obj$tau
  cEBMF.obj <- update_elbo.cEBMF(cEBMF.obj$elbo)
print(cEBMF.obj$elbo)
  return(cEBMF.obj)

}

#'@param cEBMF a cEBMF object
update_elbo.cEBMF <- function( cEBMF){
  cEBMF.obj$elbo <- c(cEBMF.obj$elbo,
                      get_objective.cEBMF (cEBMF.obj)
  )
  return(cEBMF.obj)
}


#'@param cEBMF a cEBMF object
#'@param k the factor to be updated
#'
#'
update_cEBMF <-  function(cEBMF.obj, k)
{

  #loading update
  Rk <- cal_partial_residuals.cEBMF(cEBMF.obj,k)
  l_k <- cal_expected_loading( cEBMF.obj, Rk,k)
  t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                              se      = l_k$s_i,
                              X       = cEBMF.obj$X_l )
  t_fit <- fit.mococomo(t_data)

  fitted_loading <- post_mean_sd.mococomo (t_fit )


  if( cEBMF.obj$K==1){
    cEBMF.obj$loading      <-  fitted_loading$mean
    cEBMF.obj$loading2     <- fitted_loading$sd^2+ fitted_loading$mean^2
  }else{
    cEBMF.obj$loading[,k]  <-  fitted_loading$mean
    cEBMF.obj$loading2[,k] <- fitted_loading$sd^2+ fitted_loading$mean^2
  }
  cEBMF.obj$KL_l[[k]]  <- get_KL.mococomo(t_fit)

  #factor update

  f_k <- cal_expected_factor( cEBMF.obj,Rk,k)
  t_data <- set_data_mococomo(betahat = f_k$f_j_hat,
                              se      = f_k$s_j,
                              X       = cEBMF.obj$X_f )
  t_fit <- fit.mococomo(t_data)

  fitted_factor <- post_mean_sd.mococomo (t_fit )

  if( cEBMF.obj$K==1){
    cEBMF.obj$factor      <-  fitted_factor $mean
    cEBMF.obj$factor2     <-  fitted_factor$sd^2+ fitted_factor$mean^2
  }else{
    cEBMF.obj$factor[,k]  <-  fitted_factor $mean
    cEBMF.obj$factor2[,k] <-  fitted_factor$sd^2+ fitted_factor$mean^2
  }
  cEBMF.obj$KL_f[[k]]  <- get_KL.mococomo(t_fit)

  return(cEBMF.obj)
}



#' @title prepare output of cEMBF function
#' @description prepare output of cEMBF function
#' @param cEBMF.obj a cEBMF object
#' @return a cEBMF object

out_prep.cEBMF <- function(cEBMF.obj)
{
  cEBMF.obj$Y_fit <- cal_fitted_value.cEBMF(cEBMF.obj)
}



compute_elbo.cEBMF <- function(cEBMF.obj){


  sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) + e_loglik(data, f)
}


get_objective.cEBMF = function(cEBMF.obj) {
   out <-  e_loglik  (cEBMF.obj) - sum(unlist(cEBMF.obj$KL_l)) - sum(unlist(cEBMF.obj$KL_f))

  return(out)
}


# @title Get expected log likelihood under current fit.
#
# @inheritParams flash_get_objective
#
e_loglik = function(cEBMF.obj) {
  R2    <- cal_expected_residuals.cEBMF(cEBMF.obj)
  tau   <- cEBMF.obj$tau
  e_log <-  e_loglik_R2_tau(R2, tau )
  return(e_log)
}


# adapted from flashr
e_loglik_R2_tau= function(R2, tau, data) {
  # tau can be either a scalar or a matrix:

  return(-0.5 * sum(log(2 * pi / tau) + tau * R2))
}




