
test_ibss2m<- function() {
  sim <- sim_susie(n = 100)
  fit <- ibss2m(sim$X, sim$y, L=5, tol=1e-2, estimate_prior_variance = F)
  fit2 <- binsusie(sim$X, sim$y, L=5, estimate_prior_variance = F)

  data <- with(sim, binsusie_prep_data(X, y, N, Z))
  data$shift <- 0
  data$shift_var <- 0

  fit <- data_initialize_ibss2m(data, L=5)
  fit <- update_model(fit, data)
  fit <- fit_model(fit, data)

  assignment_logp <- predict_assignment_mnsusie(fit, data)
  z <- purrr::map_int(1:nrow(assignment_logp), ~which.max(assignment_logp[.x,]))
  zsim <-  purrr::map_int(1:nrow(sim$y), ~which.max(sim$y[.x,]))

  confusion <- table(zsim, z)
  image(scale(confusion, center=F, scale=T))
  table(z, zsim)
  compute_assi
  psi <- compute_psi(fit, data)

  z <- purrr::map_int(1:nrow(data$Y), ~which.max(data$Y[.x,]))
  psi <- compute_psi(fit, data)
  pi_tilde <- sigmoid(psi)

  # logpi <- do.call(rbind, purrr::map(1:nrow(psi), ~.predict2logpi(psi[.x,])))
  # zpred <- purrr::map_int(1:nrow(logpi), ~which.max(logpi[.x,]))

  return(list(
    fit = fit,
    monotone = .monotone(fit$elbo)
  ))
}

