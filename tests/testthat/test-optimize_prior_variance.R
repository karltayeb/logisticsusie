test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test <- function(){
  n <- 1000
  betahat <- rnorm(n)
  s2 <- rep(1, n)
  lr <- 1 / exp(-0.5 * rchisq(n, 1, 0))


  opt <- optimize_prior_variance(betahat, s2, lr, pi, laplace=T, min_prior_variance=1e-5)

  f_opt <- function(x){
    purrr::map_dbl(x, ~compute_log_labf_ser(betahat, s2, lr, .x, rep(1/n, n)))
  }

  f0 <- f_opt(0)
  fopt <- f_opt(opt)
}
