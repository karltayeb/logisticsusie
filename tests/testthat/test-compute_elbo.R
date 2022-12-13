testthat::test_that("Different MoCoCoMo ELBO functions agree", {
  data <- sim_mococomo()
  fit <- fit.mococomo(data, maxiter = 10)

  K <- length(fit$f_list)
  N <- .expected_trials(fit$post_assignment)
  # update omega so jj bound is tight
  for(k in seq(K-1)){
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- N[, k]
    fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit , k=k)
    fit$logreg_list[[k]]$params$tau <- compute_tau(fit,k=k)
  }

  elbo1 <- compute_elbo.mococomo(fit)
  elbo2 <- compute_elbo2.mococomo(fit)
  elbo3 <- compute_elbo3.mococomo(fit)
  testthat::expect_equal(elbo1, elbo2)
  testthat::expect_equal(elbo1, elbo3)

})


