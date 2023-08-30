surv_uni_fun <- function(x, y, o, prior_variance, estimate_intercept = 0, ...){
  fit <- coxph(y~ x + offset(o))
  bhat <- summary(fit)$coefficients[1, 1] # bhat = -alphahat
  sd <- summary(fit)$coefficients[1, 3]
  zscore <- bhat/sd
  lr <- summary(fit)$logtest[1]/2

  lbf <- compute_lbf(zscore, sd, prior_variance)
  lbf.corr <- lbf - bhat^2/sd^2/2+ as.numeric(summary(fit)$logtest[1]/2)
  var <- compute_approx_post_var(zscore, sd, prior_variance)
  mu <- compute_approx_post_mean(var, sd, bhat)
  return(list(mu = mu, var=var, lbf=lbf.corr, intercept=0))
}
