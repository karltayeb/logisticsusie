#'@title Main FDR function
#'@description Main function gathering all the routines for different FDR estimation procedures





cFDR <- function( betahat,
                  se,
                  pvalue,
                  zscore,
                  X,
                  coverage  = 0.95,
                  maxiter   = 100,
                  tol       = 1e-3,
                  max_class = 10,
                  mult      = 2,
                  upper     = FALSE
                  )
{
  if( !(missing(betahat))&missing(se)){
    stop("Please provide standard error when providing regression estimate or use zscore input")
  }
  if(   missing(betahat)&missing(pvalue)&missing(zscore)){
    stop("Please provide one of the following entry:betahat, pvalue or zscore")
  }


  if( !missing(beta)& !missing( beta)){
    dist <- "normal"
    data <- set_data_mococomo(betahat  = betahat,
                              X = X,
                              se=  se)
    #fit mococomo model
  }
  if(!missing(pvalue)){
    dist <- "beta"
    data <- set_data_mococomo(p = pvalue,
                              X = X)
  }

  if(!missing(zscore)){
    dist <- "beta"
    data <- set_data_mococomo(p = pnorm(zscore), #using U value as in the ZAP paper
                              X = X)
  }

       fit =  fit.mococomo(data,
                           dist      = dist,
                           maxiter   = maxiter,
                           tol       = tol,
                           max_class = max_class,
                           mult      = mult,
                           upper     = upper)




       plot( est$mean, data$betahat)
return( fit)
}


prep_out_FDR_wrapper <- function(fit){
  if(dist=="normal"){
    est    <- post_mean_sd.mococomo (fit)
    lfdr   <- get_lfdr_mixnorm(fit)
    qvalue <- cal_qvalue(lfdr)
    data.frame(betahat       = fit$data$betahat,
               se            = fit$data$se,
               lfdr          = lfdr,
               qvalue        = qvalue,
               PosteriorMean = est[,1],
               PosteriorSD   = est[,2]
          )
  }

  cs <- lapply( 1:length(fit$logreg_list), function(k)
                                          get_all_cs(fit$logreg_list[[k]])
                )
}


cal_pip <-

get_lfdr_mixnorm <- function(fit){
  tt1 <-  fit$post_assignment[,1]* dnorm( fit$data$betahat, mean=0, sd= fit$data$se )
  tt2 <-    Reduce("+",lapply( 2: ncol(fit$post_assignment),
                               function(k)  fit$post_assignment[,k]*
                                            dnorm( fit$data$betahat,
                                                   mean=0,
                                                   sd= sqrt( fit$data$se^2+ fit$f_list[[k]]$var ) )))
  out <- tt1/(tt1+tt2)

  return(out)
}

cal_qvalue <- function(lfdr)
{

  torder <- order(lfdr)
  qvalue <- sapply( 1:length(lfdr), function(k )
                                    mean(lfdr[which(torder <= torder[k] )])
                    )
  return(qvalue)

}


assesor.beta <- function(fit){}
