
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
X_l <- matrix(rnorm(nrow(Y)*10), nrow= nrow(Y))
X_f <- matrix(rnorm(ncol(Y)*10), nrow= ncol(Y))


cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")




for (i in 1:5) {
  cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
}

Y_est <- Reduce("+", lapply( 1:1, function(k) cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k]) ))


testthat::test_that("cEBMF should be similar to flash without covariate", {
  expect_gt(cor(c(fitted(f)),c(Y_est)) , 0.99)
}
)

