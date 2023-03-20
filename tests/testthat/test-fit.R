###
# Multinomial susie
###

test_mn_susie <- function() {
  data <- sim_mn_susie()
  data <- with(data, mnsusie_prep_data(X, y, Z))
  data$shift <- 0
  data$shift_var <- 0
  fit <- data_initialize_sbmn_susie(data, L=5)
  fit <- update_model(fit, data, track_elbo=F)
  fit <- fit_model(fit, data)


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

