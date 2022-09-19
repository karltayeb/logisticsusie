
#rm(list = ls())
library(logisticsusie)  #Simulate data under the mococomo model
sim1  <- sim_twococomo(n=1000)


tt <- fit.mococomo(sim1)
out <- post_mean_sd.mococomo(tt)

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

plot( Y,Y_true)
X_l =data1$X

X_f =data2$X
image(Y)
image(Y_true)



K=3
u <- svd(Y,K,K)

str(u)
u <- softImpute(Y, rank.max=K,type="als", lambda=0)
plot(u$u%*%diag(u$d[1:K])%*%t(u$v),Y)
plot(u$u%*%diag(u$d[1:K])%*%t(u$v),Y_true)


dim(Y)
cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")
plot(cEBMF.obj$loading[,1],data1$betahat)
abline(a=0,b=1)
plot(cEBMF.obj$factor[,1],data2$betahat)
abline(a=0,b=1)


plot(cEBMF.obj$loading%*%t(cEBMF.obj$factor),Y)

plot(cEBMF.obj$loading%*%t(cEBMF.obj$factor),Y_true)


for ( o in 1:4){
  k=1
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


  cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )

  Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
  plot( Y_est, Y )
  abline(a=0,b=1)
}


Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
plot( Y_est, Y )
plot( Y_true, Y )
cor( c(Y_true), c(Y_est))
cor( c(Y_true), c(Y))


flash(Y)





####ash ----

library(flashr)
library(softImpute)
set.seed(1) # for reproducibility
ftrue = matrix(rnorm(200), ncol=2)
ltrue = matrix(rnorm(40), ncol=2)
ltrue[1:10, 1] = 0 # set up some sparsity
ltrue[11:20, 2] = 0
Y = ltrue %*% t(ftrue) + rnorm(2000) # set up a simulated matrix
f = flash(Y, K=1 ,ebnm_fn= 'ebnm_ash' )
ldf = f$ldf
Y_true <- ltrue %*% t(ftrue)

K=3
u <- softImpute(Y, rank.max=K,type="als", lambda=0)
plot(u$u%*%diag(u$d[1:K])%*%t(u$v),Y)
plot(u$u%*%diag(u$d[1:K])%*%t(u$v),Y_true)


dim(Y)
X_l <- matrix(rnorm(nrow(Y)*10), nrow= nrow(Y))
X_f <- matrix(rnorm(ncol(Y)*10), nrow= ncol(Y))


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

abline(a=0,b=1)


plot(cEBMF.obj$loading%*%t(cEBMF.obj$factor),Y)
abline(a=0,b=1)

plot(cEBMF.obj$loading%*%t(cEBMF.obj$factor),Y_true)
abline(a=0,b=1)

t_fit <- lm(c(Y_true)~c(Y))
summary(t_fit)

mean(f$fit$tau)
mean(cEBMF.obj$tau)
library(ashr)
for ( o in 1:4){
  k=1
  Rk <- Y
  l_k <- cal_expected_loading( cEBMF.obj, Rk,k)
  t_ash <- ash(betahat = l_k$l_i_hat,
               se      = l_k$s_i,
               mixcompdist = "normal")
  cEBMF.obj$loading[,k] <- t_ash$result$PosteriorMean
  cEBMF.obj$loading2[,k] <- t_ash$result$PosteriorSD^2 + t_ash$result$PosteriorMean^2

  f_k <- cal_expected_factor( cEBMF.obj, Rk,k)

  t_ash <- ash(betahat = f_k$f_j_hat,
               se      = f_k$s_j,
               mixcompdist = "normal")

  cEBMF.obj$factor[,k] <- t_ash$result$PosteriorMean
  cEBMF.obj$factor2[,k] <-  t_ash$result$PosteriorSD^2 + t_ash$result$PosteriorMean^2



  cEBMF.obj<- update_tau.cEBMF (cEBMF.obj )
  cEBMF.obj$tau
  Y_est <- cEBMF.obj$loading[,1]%*%t(cEBMF.obj$factor[,1])
  plot( Y_est, Y )
  abline(a=0,b=1)
}

plot( Y_est, Y_true )
f = flash(Y, K=1 ,ebnm_fn= 'ebnm_ash' )
plot( fitted(f),Y_true)
plot( Y_est, fitted(f) )


X_l <- matrix(rnorm(nrow(Y)*10), nrow= nrow(Y))
X_f <- matrix(rnorm(ncol(Y)*10), nrow= ncol(Y))


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

abline(a=0,b=1)
for ( o in 1:4){
  k=1
  Rk <- Y
  #Rk <- cal_partial_residuals.cEBMF(cEBMF.obj,k)
  l_k <- cal_expected_loading( cEBMF.obj, Rk,k)
  t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                              se      = l_k$s_i,
                              X       = cEBMF.obj$X_l )
  t_fit <- fit.mococomo(t_data)

  fitted_loading <- post_mean_sd.mococomo (t_fit )
  cEBMF.obj$loading[,k] <-  fitted_loading$mean
  cEBMF.obj$loading2[,k] <- fitted_loading$sd^2+ fitted_loading$mean^2

  #factor update

  f_k <- cal_expected_factor( cEBMF.obj,Rk,k)
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
plot( fitted(f),Y_true)
plot( Y_est, fitted(f) )


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")



####Hand run update----

Y_est <- cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k])
plot( Y_est, Y_true )

cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

for (i in 1:5) {
  cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
}
Y_est <- Reduce("+", lapply( 1:cEBMF.obj$K, function(k) cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k]) ))


plot( Y_est, fitted(f) )
abline(a=0,b=1)
plot( Y_est, Y_true)
points(fitted(f),Y_true  , col="green" )
