set.seed(1)
library(logisticsusie)
sim  <- sim_mococomo_beta(n=100)
#preparing the data
data <- set_data_mococomo(p = sim$p,
                          X = sim$X)
#fit mococomo model
maxiter   = 100
tol       = 1e-3
max_class = 10
mult      = 2
upper     = TRUE


#working init
fit <- init.mococomo(data,
                     max_class = max_class,
                     mult   = mult,
                     upper     = upper
)
#working elbo
testthat::test_that("agreement between ELBO",{
  testthat::expect_equal(compute_elbo.mococomo(fit),compute_elbo2.mococomo(fit))
  testthat::expect_equal(compute_elbo.mococomo(fit),compute_elbo3.mococomo(fit))
  }
)

#this should work
#works with two distribution (left and right)
tfit <- fit.mococomo  (data,upper=TRUE)
#works with a singel left distribution
tfit <- fit.mococomo  (data,upper=FALSE)
