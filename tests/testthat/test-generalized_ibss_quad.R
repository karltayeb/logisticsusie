test_that("multiplication works", {
  sim <- sim_ser()
  fit0 <- with(sim, generalized_ibss(X, y, L=5))
  fit0 <- with(sim, generalized_ibss(X, y, L=1))
  fit1 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F))
  fit2 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F, n=1024))
  plot(fit$fits[[1]]$lbf, fit0$fits[[1]]$lbf); abline(0, 1, col='red')
})


test_that("estimate prior variance works", {
  sim <- sim_ser()
  fiteb <- with(sim, generalized_ibss(X, y, L=1, estimate_prior_variance=T))
  fit <- with(sim, generalized_ibss(X, y, L=1, estimate_prior_variance=F))

  fit0 <- with(sim, generalized_ibss(X, y, L=1))
  fit1 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F))
  fit2 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F, n=1024))
  plot(fit$fits[[1]]$lbf, fit0$fits[[1]]$lbf); abline(0, 1, col='red')
})
