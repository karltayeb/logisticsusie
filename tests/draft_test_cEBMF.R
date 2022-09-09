library(flashr)
#rm(list = ls())
library(logisticsusie)  #Simulate data under the mococomo model
sim1  <- sim_twococomo(n=1000)


tt <- fit.mococomo(sim1)
out <- post_mean_sd.mococomo(fit)

plot( out$mean, sim1$beta)
#preparing the data
data1 <-set_data_mococomo(betahat = sim1$betahat,
                          se = sim1$se ,
                          X  = sim1$X)
sim2  <- sim_twococomo(20)
#preparing the data
data2 <-set_data_mococomo(betahat = sim2$betahat,
                          se = sim2$se ,
                          X  = sim2$X)



Y_true <- sim1$beta%*%t(sim2$beta)
Y <- data1$betahat%*%t(data2$betahat)
X_l =data1$X

X_f =data2$X
image(Y)
image(Y_true)


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=1)




for ( o in 1:4){
  k=1
  l_k <- cal_expected_loading( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                              se      = l_k$s_i,
                              X       = cEBMF.obj$X_l )
  t_fit <- fit.mococomo(data =t_data)
  str(t_fit)

  l_i_fit <-  post_mean_sd.mococomo (t_fit )$mean
  plot(l_i_fit,data1$betahat)


  cEBMF.obj$loading[,k] <- l_i_fit



  f_k <- cal_expected_factor( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = f_k$f_j_hat,
                              se      = f_k$s_j,
                              X       = cEBMF.obj$X_f )
  t_fit <- fit.mococomo(t_data)
  f_j_fit <- post_mean_sd.mococomo (t_fit )$mean
  plot(f_j_fit,data2$betahat)

  cEBMF.obj$loading[,k] <-  post_mean_sd.mococomo (fit )$mean

  cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )



}


Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
plot( Y_est, Y )
plot( Y_true, Y )
cor( c(Y_true), c(Y_est))
cor( c(Y_true), c(Y))






#### tesst using ash

library(flashr)
#rm(list = ls())
library(logisticsusie)  #Simulate data under the mococomo model
sim1  <- sim_twococomo(n=1000)


tt <- fit.mococomo(sim1)
out <- post_mean_sd.mococomo(fit)

plot( out$mean, sim1$beta)
#preparing the data
data1 <-set_data_mococomo(betahat = sim1$betahat,
                          se = sim1$se ,
                          X  = sim1$X)
sim2  <- sim_twococomo(20)
#preparing the data
data2 <-set_data_mococomo(betahat = sim2$betahat,
                          se = sim2$se ,
                          X  = sim2$X)



Y_true <- sim1$beta%*%t(sim2$beta)
Y <- data1$betahat%*%t(data2$betahat)
X_l =data1$X

X_f =data2$X
image(Y)
image(Y_true)


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=1)




for ( o in 1:4){
  k=1
  l_k <- cal_expected_loading( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                              se      = l_k$s_i,
                              X       = cEBMF.obj$X_l )

  t_fit <- ash(t_data$betahat,t_data$se)
  l_i_fit <-t_fit$result$PosteriorMean

  plot(l_i_fit,data1$betahat)


  cEBMF.obj$loading[,k] <- l_i_fit



  f_k <- cal_expected_factor( cEBMF.obj,k)
  t_data <- set_data_mococomo(betahat = f_k$f_j_hat,
                              se      = f_k$s_j,
                              X       = cEBMF.obj$X_f )
  t_fit <- ash(t_data$betahat,t_data$se)
  f_j_fit <- t_fit$result$PosteriorMean
  plot(f_j_fit,data2$betahat)

  cEBMF.obj$loading[,k] <- f_j_fit

  cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )



}


Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
plot( Y_est, Y )
plot( Y_true, Y )
cor( c(Y_true), c(Y_est))
cor( c(Y_true), c(Y))

