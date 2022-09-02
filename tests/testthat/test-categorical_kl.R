test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


.generate_pi <- function(K=50){
  logits <- rnorm(K)
  logits <- logits - logSumExp(logits)
  pi <- exp(logits)
  assertthat::are_equal(sum(pi), 1)
  return(pi)
}

testthat::test_that("categorical KL is positive", {
  set.seed(1)
  for (i in seq(100)){
    pi1 <- .generate_pi()
    pi2 <- .generate_pi()
    testthat::expect_true(categorical_kl(pi1, pi2) >=0)
  }
})
