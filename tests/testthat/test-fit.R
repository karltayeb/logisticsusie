###
# Multinomial susie
###

test_mn_susie <- function() {
  data <- sim_mn_susie()
  fit <- fit.mnsusie(data, L = 3, maxiter = 50)
  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}


testthat::test_that("Multinomial SuSiE", {
  for (i in seq(5)) {
    testthat::expect_true(test_mn_susie()$monotone)
  }
})

