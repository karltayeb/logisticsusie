###
# Multinomial susie
###

test_mn_susie <- function() {
  sim <- sim_mn_susie()
  data <- with(sim, mnsusie_prep_data(X, y, Z))
  data$shift <- 0
  data$shift_var <- 0
  fit <- data_initialize_sbmn_susie(data, L=5)
  fit <- update_model(fit, data, track_elbo=F)
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

testthat::test_that("Multinomial SuSiE", {
  for (i in seq(5)) {
    testthat::expect_true(test_mn_susie()$monotone)
  }
})


test_cb_susie <- function() {
  sim <- sim_mn_susie()
  data <- with(sim, mnsusie_prep_data(X, y, Z))
  data$shift <- 0
  data$shift_var <- 0

  # cb
  fit <- data_initialize_cb_susie(data, L=5)
  fit <- update_model(fit, data, track_elbo=F)
  fit <- fit_model(fit, data)

  assignment_logp <- predict_log_assignment_cbsusie(fit, data)
  z <- purrr::map_int(1:nrow(assignment_logp), ~which.max(assignment_logp[.x,]))
  zsim <-  purrr::map_int(1:nrow(sim$y), ~which.max(sim$y[.x,]))
  confusion <- table(zsim, z)
  image(scale(confusion, center=F, scale=T))

  # sb
  fit2 <- data_initialize_sbmn_susie(data, L=5)
  fit2 <- update_model(fit2, data, track_elbo=F)
  fit2 <- fit_model(fit2, data)

  assignment_logp2 <- predict_log_assignment_sbnmsusie(fit2, data)
  z2 <- purrr::map_int(1:nrow(assignment_logp2), ~which.max(assignment_logp2[.x,]))


  zsim <-  purrr::map_int(1:nrow(sim$y), ~which.max(sim$y[.x,]))
  confusion <- table(zsim, z)
  confusion2 <- table(zsim, z2)

  image(scale(confusion, center=F, scale=T))
  image(scale(confusion2, center=F, scale=T))

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
