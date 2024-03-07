make_c1_sim <- function(){
  ##### Minimal working example
  c1 <- gseasusie::load_msigdb_geneset_x('C1')

  # sample random 5k background genes
  set.seed(0)
  background <- sample(rownames(c1$X), 5000)

  # sample GSs to enrich, picked b0 so list is ~ 500 genes
  enriched_gs <- sample(colnames(c1$X), 3)
  b0 <- -2.2
  b <- 3 *abs(rnorm(length(enriched_gs)))
  logit <- b0 + (c1$X[background, enriched_gs] %*% b)[,1]
  y1 <- rbinom(length(logit), 1, 1/(1 + exp(-logit)))
  list1 <- background[y1==1]
  X <- c1$X[background,]
  X <- X[, Matrix::colSums(X) > 1]
  idx <- which(colnames(X) %in% enriched_gs)

  list(X = X, y = y1, idx = idx, b = b, b0 = b0, logits=logit)
}


test_that("c1 simulation, slow but accurate", {
  sim <- make_c1_sim()
  serfit <- with(sim, logistic_ser(X, y, prior_variance = 1)) # ~30 seconds
  gibssfit <- with(sim, ibss_from_ser(X, y, L=5, maxit=10, ser_function = logistic_ser)) #~15 minutes, 6 iterations
  saveRDS(gibssfit, 'c1_sim_gibssfit.rds')
})


test_that("c1 simulation, slow but accurate", {
  set.seed(0)
  sim <- sim_ser(n=50000)

  x <- sim$X[,1]
  y <- sim$y

  # what we currently do-- center quadrature on
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=1)$logZ # -17305.19
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=5)$logZ # -17305.19
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=100)$logZ # -17711.58 WRONG
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=Inf)$logZ # -17305.19 CORRECT when we use infinite interval

  # for very peaked functions, `integrate` can have trouble finding the non-zero region
  logistic_bayes(x, y, o = rep(0, length(y)), eps = -1)$logZ # -17381.66 WRONG
  logistic_bayes(x, y, o = rep(0, length(y)), eps = -1, width = Inf)$logZ # -17305.19 CORRECT when we use infinite interval

  sim <- make_c1_sim()

  x <- sim$X[,66]
  y <- sim$y

  # what we currently do-- center quadrature on
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=1)$logZ # -1660.454 WRONG
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=5)$logZ # -1660.391
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=100)$logZ # -1660.391 CORRECT
  logistic_bayes(x, y, o = rep(0, length(y)), eps = 0, width=Inf)$logZ # -1660.391 CORRECT

  # for very peaked functions, `integrate` can have trouble finding the non-zero region
  logistic_bayes(x, y, o = rep(0, length(y)), eps = -3.506314)$logZ # -1660.397 off a little
  logistic_bayes(x, y, o = rep(0, length(y)), eps = -3.506314, width = Inf)$logZ # -1660.391 CORRECT when we use infinite interval

  ridgefit <- ridge(x, y, prior_variance=66)
  ridge_intercept <- with(ridgefit, list(intercept = intercept, tau0 = tau0))
  f <- make_log_joint_ridge(x, y, rep(0, length(y)), ridge_intercept)
  g <- make_laplace_approx(f, ridgefit)
  C <- f(3.450823) # optimize f

  f2 <- Exp(Shift(f, C))
  f3 <- Exp(f)
  g2 <- Exp(Shift(g, C))

  ridgefit$logZ # Laplace
  log(integrate(g2, -Inf, Inf)$value) + C # confirm Laplace
  log(integrate(f2, -Inf, Inf)$value) + C # CORRECT

  log(integrate(f3, -Inf, Inf)$value) # CORRECT
  C  + 0.5 * log(2 * pi / ridgefit$tau)
  ridgefit$logZ

  log(integrate(f2, -10, 10)$value) + C # CORRECT
  log(integrate(f2, -100, 100)$value) + C # CORRECT
  log(integrate(f2, -1000, 1000)$value) + C # CORRECT


  # now we don't give it the optimum value
  ridgefit <- ridge(x, y)
  ridge_intercept <- with(ridgefit, list(intercept = intercept, tau0 = tau0))
  f <- make_log_joint_ridge(x, y, rep(0, length(y)), ridge_intercept)

  f3 <- Exp(f)
  log(integrate(f3, -10, 10)$value) # -Inf
  log(integrate(f3, -100, 100)$value) # NAN
  log(integrate(f3, -1000, 1000)$value) # NAN
  log(integrate(f3, -Inf, Inf)$value) # NAN

  sim <- make_c1_sim()
  serfit <- with(sim, logistic_ser(X, y, prior_variance = 1)) # ~30 seconds
  gibssfit <- with(sim, ibss_from_ser(X, y, L=5, maxit=10, ser_function = logistic_ser)) #~15 minutes, 6 iterations
  saveRDS(gibssfit, 'c1_sim_gibssfit.rds')
})

test_that("2d", {
  sim <- make_c1_sim()
  x <- sim$X[,66]
  y <- sim$y

  f <- make_log_joint2d(x, y, rep(0, length(y)), 1)
  b <- matrix(c(-2.15, 3.7), ncol=1)
  f(b)

  C <- f(b)[1, 1]
  f2 <- function(b){ exp(f(b) - C) }

  cubature::pcubature(
    f = f2,
    lowerLimit = rep(-10, 2),
    upperLimit = rep(10, 2),
    vectorInterface = T
  )

  f(matrix(c(-2.15, 3.7), ncol=1))
})


