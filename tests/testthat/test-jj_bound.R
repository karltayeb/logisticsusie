testthat::test_that("binomial and logistic bound agree", {
  data <- sim_ser()
  fit <- fit.binser(data)
  testthat::expect_equal(jj_bound.binser(fit), jj_bound.logistic(fit))
})

testthat::test_that("Explicit ELBO and  JJ bound agree", {
  data <- sim_ser()
  fit <- fit.binser(data)

  a <- jj_bound.binser(fit)
  b <- explicit_elbo.binser(fit)
  testthat::expect_equal(a, b)
})


testthat::test_that("MOCOMO JJ = sum(Bin SuSiE JJ)", {
  data <- sim_mococomo()
  fit <- fit.mococomo(data, maxiter = 10)

  K <- length(fit$f_list)
  N <- .expected_trials(fit$post_assignment)

  # update omega so jj bound is tight
  for(k in seq(K-1)){
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- N[, k]
    fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit$logreg_list[[k]])
    fit$logreg_list[[k]]$params$tau <- compute_tau(fit$logreg_list[[k]])
  }

  jj1 <- rowSums(compute_assignment_jj_bound.mococomo(fit) * fit$post_assignment)
  jj2 <- rowSums(do.call(cbind, purrr::map(fit$logreg_list, jj_bound.binsusie)))
  testthat::expect_equal(jj1, jj2)
})
