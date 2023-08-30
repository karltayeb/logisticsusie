surv_uni_fun <- function(x, y, o, prior_variance, estimate_intercept = 0, ...){
  fit <- survival::coxph(y ~ x + offset(o))
  smry <- summary(fit)
  intercept <- 0
  bhat <- smry$coefficients[1, 1] # bhat = -alphahat
  sd <- smry$coefficients[1, 3]
  lr <- smry$logtest[1]/2

  res <- list(betahat = bhat, shat2 = sd**2, intercept = 0, lr = lr)
  return(res)
}

survival_ser <- asymptotic_ser_from_univariate(surv_uni_fun)

survival_ibss <- function(X, y, L, ...){
  res <- ibss_from_ser(X, y, L, ser_fun = survival_ser, ...)
}
