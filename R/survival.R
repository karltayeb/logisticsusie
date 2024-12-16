surv_uni_fun <- function (x, y, e, prior_variance,
                          estimate_intercept = 0, ...) {
  v0   <- prior_variance						 
  fit  <- coxph(y ~ x + offset(e))
  out  <- summary(fit)$coefficients
  bhat <- out[1,"coef"]
  s    <- out[1,"se(coef)"]
  z    <- bhat/s
  lbf  <- log(s^2/(v0 + s^2))/2 + z^2/2*v0/(v0 + s^2)
  lbf  <- lbf - z^2/2 + summary(fit)$logtest["test"]/2
  v1   <- 1/(1/v0 + 1/s^2)
  mu1  <- v1*bhat/s^2
  return(list(mu = mu1,var = v1,lbf = lbf,intercept = 0))
}

survival_ibss <- function(X, y, L, ...){
  res <- ibss_from_ser(X, y, L, ser_fun = surv_uni_fun, ...)
}
