
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


t1 <- init_cEBMF (Y, X_l,X_f,K=1, init_type = "udv_si")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t2 <- init_cEBMF (Y, X_l,X_f,K=5, init_type = "udv_si")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t3 <- init_cEBMF (Y, X_l,X_f,K=15, init_type = "udv_si")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t1 <- init_cEBMF (Y, X_l,X_f,K=1, init_type = "udv")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t2 <- init_cEBMF (Y, X_l,X_f,K=5, init_type = "udv")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t3 <- init_cEBMF (Y, X_l,X_f,K=15, init_type = "udv_si_svd")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t1 <- init_cEBMF (Y, X_l,X_f,K=1, init_type = "udv_si_svd")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t2 <- init_cEBMF (Y, X_l,X_f,K=5, init_type = "udv_si_svd")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)

t3 <- init_cEBMF (Y, X_l,X_f,K=15, init_type = "udv_si_svd")
cEBMF.obj <- cEBMF_iter(cEBMF.obj)




Y_est <- cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k])
plot( Y_est, Y_true )

cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")

for (i in 1:5) {
  cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
  print(i)
}
Y_est <- Reduce("+", lapply( 1:cEBMF.obj$K, function(k) cEBMF.obj$loading[,k]%*%t(cEBMF.obj$factor[,k]) ))


plot(cEBMF.obj$elbo)
testthat::test_that("exemple elbo should be monotonic",{
  expect_equal(.monotone(cEBMF.obj$elbo), TRUE)
}
)

testthat::test_that("exemple elbo is not constant",{
  expect_gt(var(cEBMF.obj$elbo),0)
}
)
cEBMF.obj <- init_cEBMF (Y, X_l,X_f,K=3, init_type = "udv_si")
cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
l1 <- cEBMF.obj$loading[,1]
f1 <- cEBMF.obj$loading[,1]
cEBMF.obj <- cEBMF_iter  (cEBMF.obj)
l2 <- cEBMF.obj$loading[,1]
f2 <- cEBMF.obj$loading[,1]


plot( l1, l2)
(sum(!(l1==l2)))


testthat::test_that("the loading are updated",{
  expect_gt((sum(!(l1==l2))),0)
}
)

testthat::test_that("the loading are updated",{
  expect_gt((sum(!(f1==f2))),0)
}
)
