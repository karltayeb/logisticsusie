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


testthat::test_that('Mococomo bound can be computed from SuSiE bounds', {
  # confirm that
  # \sum_k(JJ + KL(omega)) = E[logp(y | \beta, \omega)]
  data <- sim_mococomo()
  fit <- fit.mococomo(data, maxiter = 5)

  K <- length(fit$f_list)
  N <- .expected_trials(fit$post_assignment)
  # update omega so jj bound is tight
  for(k in seq(K-1)){
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- N[, k]
    fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit$logreg_list[[k]])
    fit$logreg_list[[k]]$params$tau <- compute_tau(fit$logreg_list[[k]])
  }


  logp <- compute_assignment_loglikelihood.mococomo(fit)
  post_assignment <- fit$post_assignment
  a <- rowSums(logp * post_assignment)

  # compute E[logp(z | beta, omega)] as sum of jj bounds
  jj <- do.call(cbind, purrr::map(fit$logreg_list, jj_bound.binsusie))
  omega_kl <- do.call(cbind, purrr::map(fit$logreg_list, function(x) pg_kl(x$data$N, x$params$xi)))
  b <- rowSums(jj + omega_kl)

  testthat::expect_equal(a, b)
})


testthat::test_that('Mococomo bound can be computed from SuSiE bounds 2', {
  # confirm that
  # \sum_k(JJ + KL(omega)) = E[logp(y | \beta, \omega)]
  data <- sim_mococomo()
  fit <- fit.mococomo(data, maxiter = 5)

  K <- length(fit$f_list)
  N <- .expected_trials(fit$post_assignment)
  # update omega so jj bound is tight
  for(k in seq(K-1)){
    fit$logreg_list[[k]]$data$y <- fit$post_assignment[, k]
    fit$logreg_list[[k]]$data$N <- N[, k]
    fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit$logreg_list[[k]])
    fit$logreg_list[[k]]$params$tau <- compute_tau(fit$logreg_list[[k]])
  }


  logp <- compute_assignment_jj_bound.mococomo(fit)
  post_assignment <- fit$post_assignment
  a <- rowSums(logp * post_assignment)

  # compute E[logp(z | beta, omega)] as sum of jj bounds
  jj <- do.call(cbind, purrr::map(fit$logreg_list, jj_bound.binsusie))
  b <- rowSums(jj)

  testthat::expect_equal(a, b)
})
