library(flashr)
#rm(list = ls())
library(logisticsusie)  #Simulate data under the mococomo model
sim  <- sim_twococomo(n=1000)
#preparing the data
data1 <-set_data_mococomo(betahat = sim$betahat,
                          se = sim$se ,
                          X  = sim$X)
sim  <- sim_twococomo(20)
#preparing the data
data2 <-set_data_mococomo(betahat = sim$betahat,
                          se = sim$se ,
                          X  = sim$X)

Y <- data1$betahat%*%t(data2$betahat)
X_l =data1$X

X_f =data2$X
image(Y)
cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=1)
k=1
l_k <- cal_expected_loading( cEBMF.obj,k)
t_data <- set_data_mococomo(betahat = l_k$l_i_hat,
                            se      = l_k$s_i,
                            X       = cEBMF.obj$X_l )
t_fit <- fit.mococomo(data =t_data)
str(t_fit)

l_i_fit <- post_mean.mococomo (fit , t_data)
t_fit$logreg_list$var

update_loading <- function(cEBMF.obj, k,  mococomo.obj)
{
  cEBMF.obj$loading[,k] <-  t_fit$logreg_list
}
