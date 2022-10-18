
#rm(list = ls())
library(logisticsusie)  #Simulate data under the mococomo model
sim11  <- sim_ordinal_mod(n=10000,beta=0,se=2,alpha_start = 0)


#preparing the data

data11 <-set_data_mococomo(betahat = c(rep(0,10),sim11$obs),
                          se = rep(sim11$se,(length( sim11$obs)+10)) ,
                          X  =rbind(0*sim11$X[1:10,],sim11$X))

res <- fit.mococomo(data11)
res$post_assignment
library(ashr)
res_ash <- ash(data11$betahat, data11$se, mixcompdist="normal")

t_out <- post_mean_sd.mococomo(res)

plot( sim11$true_obs, t_out[-c(1:10),1])
points( res_ash$result$PosteriorMean, t_out[,1], col="green")
abline(a=0,b=1)


sim12  <- sim_twococomo(n=100, beta=2)

fit <- fit.mococomo(data11)
data12 <-set_data_mococomo(betahat = sim12$betahat,
                           se = sim12$se ,
                           X  = sim12$X)
sim21  <- sim_twococomo(20)
#preparing the data
data21 <-set_data_mococomo(betahat = sim21$betahat,
                          se = sim21$se ,
                          X  = sim21$X)
sim22  <- sim_twococomo(20)
data22 <-set_data_mococomo(betahat = sim22$betahat,
                          se = sim22$se ,
                          X  = sim22$X)

Y_true <- sim11$beta%*%t(sim21$beta)+ sim12$beta%*%t(sim22$beta)
Y <- data11$betahat%*%t(data21$betahat) +data12$betahat%*%t(data22$betahat)

plot( Y,Y_true)
X_l =data11$X

X_f =data22$X
image(Y)
image(Y_true)



K=3
dim(Y)
library(softImpute)
cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")
plot(cEBMF.obj$loading[,1] ,data11$betahat)



for ( o in 1:10){
  for ( k in 1:K){
    Rk <- cal_partial_residuals.cEBMF(cEBMF.obj,k)
    l_k <- cal_expected_loading( cEBMF.obj, Rk,k)
    t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                                se      = l_k$s_i,
                                X       = cEBMF.obj$X_l )
    t_fit <- fit.mococomo(t_data)

    fitted_loading <- post_mean_sd.mococomo (t_fit )
    cEBMF.obj$loading[,k] <-  fitted_loading$mean
    cEBMF.obj$loading2[,k] <- fitted_loading$sd^2+ fitted_loading$mean^2

    #factor update

    f_k <- cal_expected_factor( cEBMF.obj, Rk,k)
    t_data <- set_data_mococomo(betahat = f_k$f_j_hat,
                                se      = f_k$s_j,
                                X       = cEBMF.obj$X_f )
    t_fit <- fit.mococomo(t_data)

    fitted_factor <- post_mean_sd.mococomo (t_fit )
    cEBMF.obj$factor[,k] <-  fitted_factor $mean
    cEBMF.obj$factor2[,k] <-  fitted_factor$sd^2+ fitted_factor$mean^2


    cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )

    Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
    plot( Y_est, Y )
    abline(a=0,b=1)
  }

}


Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])+cEBMF.obj$loading[,2]%*%t(cEBMF.obj$factor[,2])
library(flashr)
f <- flash(Y)
plot( Y_est, Y )
plot( Y_true, Y )
plot( Y_true, Y_est )
points(Y_true, fitted(f), col="green")
cor( c(Y_true), c(Y_est))
cor( c(Y_true), c(Y))
cor( c(Y_true), c(fitted(f)))










####Hand run update----

library(flashr)
library(softImpute)
ftrue = matrix(rnorm(200), ncol=2)
ltrue = matrix(rnorm(40), ncol=2)
ltrue[1:10, 1] = 0 # set up some sparsity
ltrue[11:20, 2] = 0
Y = ltrue %*% t(ftrue) + rnorm(2000) # set up a simulated matrix
f = flash(Y, K=3 ,#ebnm_fn= 'ebnm_ash' ,
          var_type = "constant")
ldf = f$ldf
Y_true <- ltrue %*% t(ftrue)
X_l <- matrix(rnorm(nrow(Y)*10), nrow= nrow(Y))
X_f <- matrix(rnorm(ncol(Y)*10), nrow= ncol(Y))

cEBMF.obj <- cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

cEBMF.obj$elbo


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

for (i in 1:10) {
  cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
 # print(i)

  print(cbind(c(cEBMF.obj$KL_f),c(cEBMF.obj$KL_l)))
  print(mean(cEBMF.obj$tau))
  print(cEBMF.obj$elbo)

}



Y_est <- Reduce("+", lapply( 1:cEBMF.obj$K, function(k) cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k]) ))
f = flash(Y, K=3 ,#ebnm_fn= 'ebnm_ash' ,
          var_type = "constant")
plot( cEBMF.obj$elbo)
plot( Y_est,Y)
 points( Y_est,Y_true, col="green")

 plot( fitted(f),Y_true)
 points( Y_est,Y_true, col="green")
 plot( Y_est,fitted(f))


 cEBMF.fit <- cEBMF (Y, X_l,X_f,K=1, init_type = "udv_si")
 cEBMF.fit$elbo

 cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=1, init_type = "udv_si")

 for (i in 1:3) {
   cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
   # print(i)

   print(cbind(c(cEBMF.obj$KL_f),c(cEBMF.obj$KL_l)))
   print(mean(cEBMF.obj$tau))
   print(cEBMF.obj$elbo)

 }
 Y_est <- Reduce("+", lapply( 1:cEBMF.obj$K, function(k) cEBMF.obj$loading %*%t(cEBMF.obj$factor ) ))
 Y_fit <- Reduce("+", lapply( 1:cEBMF.fit $K, function(k) cEBMF.fit $loading %*%t(cEBMF.fit $factor ) ))
 plot( Y_est, cEBMF.fit$Y_fit)
 abline(a=0,b=1)
 f = flash(Y, K=1 ,#ebnm_fn= 'ebnm_ash' ,
           var_type = "constant")
 plot( cEBMF.obj$elbo)
 plot( Y_est,Y)
 points( Y_est,Y_true, col="green")

 plot( fitted(f),Y_true)
 points( Y_est,Y_true, col="green")
 plot( Y_est,fitted(f))
