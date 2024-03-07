test_that("multiplication works", {
  sim <- sim_ser()
  fit0 <- with(sim, generalized_ibss(X, y, L=5))

  psi <- (sim$X %*% (fit0$alpha * fit0$mu)[1,] )[,1]
  with(sim, fit_glm_ser(X, y, o=psi))

  glm(sim$y ~ 1 + offset(psi), family='binomial')

  fit0 <- with(sim, generalized_ibss(X, y, L=1))
  fit1 <- with(sim, generalized_ibss_quad(X, y, L=5, maxit=10, verbose=F))
  fit2 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F, n=1024))
  plot(fit$fits[[1]]$lbf, fit0$fits[[1]]$lbf); abline(0, 1, col='red')
})

test_that("quad ser", {
  sim <- sim_ser()
  fit0 <- with(sim, generalized_ibss(X, y, L=5))

  quadser <- with(sim, fit_quad_ser(X, y))

  with(sim, ibss_from_ser(X, y, L=5, ser_function = fit_quad_ser, maxit=10))

  psi <- (sim$X %*% (fit0$alpha * fit0$mu)[1,] )[,1]
  with(sim, fit_glm_ser(X, y, o=psi))

  glm(sim$y ~ 1 + offset(psi), family='binomial')

  fit0 <- with(sim, generalized_ibss(X, y, L=1))
  fit1 <- with(sim, generalized_ibss_quad(X, y, L=5, maxit=10, verbose=F))
  fit2 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F, n=1024))
  plot(fit$fits[[1]]$lbf, fit0$fits[[1]]$lbf); abline(0, 1, col='red')
})


test_that("estimate prior variance works", {
  sim <- sim_ser()
  fiteb <- with(sim, generalized_ibss(X, y, L=1, estimate_prior_variance=T))
  fit <- with(sim, generalized_ibss(X, y, L=1, estimate_prior_variance=F))

  fit0 <- with(sim, generalized_ibss(X, y, L=1))
  fit1 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F))
  fit2 <- with(sim, generalized_ibss_quad(X, y, L=1, maxit=1, verbose=F, n=1024))
  plot(fit$fits[[1]]$lbf, fit0$fits[[1]]$lbf); abline(0, 1, col='red')
})

test_that("GSEA example", {
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


  # gene set matrix restricted to background, keep non-empty gene sets
  X <- c1$X[background,]
  X <- X[, Matrix::colSums(X) > 1]

  serfit1 <- fit_glm_ser(X, y1)

  augmented_data <- augment_binary_data(X, y1, a = 1.)
  serfit2 <- with(augmented_data, fit_glm_ser(X, y))

  fit <- with(augmented_data, generalized_ibss(X, y, L=5, estimate_prior_variance=T, maxit=10))

  fit$prior_variance

  p <- length(serfit$betahat)
  pi <- rep(1/p, p)

  with(serfit, max(compute_log_labf(betahat, shat2, lr, 19.)))

  serfit <- fit_glm_ser(X, y1, estimate_prior_variance = T)
  labf <- with(serfit, compute_log_labf(betahat, shat2, lr, 1.))
  z <- serfit$betahat/sqrt(serfit$shat2)^2
  plot(labf, z)
  with(serfit, plot(lr, z))

  with(serfit, compute_log_labf_ser(betahat, shat2, lr, prior_variance=10., pi = pi))

  f_opt <- function(prior_variance){
    with(serfit, {
      p <- length(betahat)
      compute_log_labf_ser(betahat, shat2, lr, prior_variance, rep(1/p, p))
    })
  }
  f <- Vectorize(f_opt)
  plot(f, xlim = c(0, 50))

  par(mfrow=c(1, 2))

  serfit <- fit_glm_ser(X, y1, estimate_prior_variance = T)
  labf <- with(serfit, compute_log_labf(betahat, shat2, lr, 1.))
  z <- serfit$betahat/sqrt(serfit$shat2)^2
  bad_points <- abs(z < 0.001) & (serfit$lr > 1)
  good_points <- !bad_points
  plot(labf[good_points], z[good_points], xlim=range(labf), ylim=range(z),
       xlab='logBF (Laplace Approximation)', ylab = 'z-score')
  points(labf[bad_points], z[bad_points], col='red')


  augmented_data <- augment_binary_data(X, y1, 1.)
  serfitaug <- with(augmented_data, fit_glm_ser(X, y, estimate_prior_variance=T))
  labf <- with(serfitaug, compute_log_labf(betahat, shat2, lr, 1.))
  z <- with(serfitaug, betahat/sqrt(shat2)^2)
  plot(labf[good_points], z[good_points], xlim=range(labf), ylim=range(z),
       xlab='logBF (Laplace Approximation)', ylab = 'z-score')
  points(labf[bad_points], z[bad_points], col='red')

  augmented_data <- augment_binary_data(X, y1, 10.)
  serfitaug <- with(augmented_data, fit_glm_ser(X, y, estimate_prior_variance=T))
  labf <- with(serfitaug, compute_log_labf(betahat, shat2, lr, 1.))
  z <- with(serfitaug, betahat/sqrt(shat2)^2)
  plot(labf, z)

})


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

test_that("C1 SER comparison", {
  library(tictoc)
  sim <- make_c1_sim()
  # estimate prior variance with UVB
  f <- function(prior_variance){
    with(sim, fit_uvb_ser(X, y, prior_variance=prior_variance, estimate_prior_variance=F))$lbf_model
  }
  #prior_variance <- optimize(f, c(1e-10, 100), maximum = T)$maximum
  prior_variance = 66.63581
  #plot(prior_variance*seq(0.5, 2, 0.25), Vectorize(f)(prior_variance*seq(0.5, 2, 0.25)))

  tic()
  uvbser <- with(sim, fit_uvb_ser(X, y, prior_variance=prior_variance, estimate_prior_variance=F))
  toc()

  # estimate prior variance with Laplace MLE
  glmser <- with(sim, fit_glm_ser(X, y, prior_variance=1, estimate_prior_variance=T))
  abfser <- with(sim, fit_glm_ser(X, y, prior_variance=1, estimate_prior_variance=T, laplace = F))

  glue::glue('UVB estimated prior variance = {uvbser$prior_variance}\n Laplace estimated prior variance = {glmser$prior_variance}\n Wakefield estimated prior variance = {abfser$prior_variance}')

  # LAPLACE
  f2 <- function(var0){with(glmser, compute_log_labf(betahat, shat2, lr, var0))}
  sigma2_grid <- seq(0.1, 50, by=0.1)
  laplace_log_bfs <- do.call('cbind', purrr::map(sigma2_grid, f2))

  f3 <- function(var0){with(glmser, compute_log_labf_ser(betahat, shat2, lr, var0, pi = 1/length(betahat)))}
  laplace_log_bf_sers <- purrr::map_dbl(sigma2_grid, f3)

  # here we plot the in the laplace approximait
  par(mfrow=c(1,2))
  plot(sigma2_grid, laplace_log_bfs[1,], ylim = c(-1, 1) + range(laplace_log_bfs), type = 'l',  ylab = 'Laplace log(BF)', xlab = 'prior variance', main='Individual variables')
  for(i in 2:nrow(laplace_log_bfs)){
    lines(sigma2_grid, laplace_log_bfs[i,], ylim = c(-1, 1) + range(laplace_log_bfs), type = 'l')
  }
  lines(sigma2_grid, laplace_log_bf_sers, col='red')
  plot(sigma2_grid, laplace_log_bf_sers, col='red', type='l', ylab = 'Laplace log(BF_SER)', xlab = 'prior variance', main='SER')


  # ABF
  f4 <- function(var0){with(glmser, compute_log_abf(betahat, shat2, var0))}
  sigma2_grid <- seq(0.1, 50, by=0.1)
  log_abfs <- do.call('cbind', purrr::map(sigma2_grid, f4))

  f5 <- function(var0){with(glmser, compute_log_abf_ser(betahat, shat2, var0, pi = 1/length(betahat)))}
  log_abf_ser <- purrr::map_dbl(sigma2_grid, f5)

  # here we plot the in the laplace approximait
  par(mfrow=c(1,2))
  plot(sigma2_grid, laplace_log_bfs[1,], ylim = c(-1, 1) + range(log_abfs), type = 'l',  ylab = 'Wakefield log(BF)', xlab = 'prior variance', main='Individual variables')
  for(i in 2:nrow(log_abfs)){
    lines(sigma2_grid, log_abfs[i,], ylim = c(-1, 1) + range(log_abfs), type = 'l')
  }
  lines(sigma2_grid, log_abf_ser, col='red')
  plot(sigma2_grid, log_abf_ser, col='red', type='l', ylab = 'Wakefield log(BF_SER)', xlab = 'prior variance', main='SER')


  sigma2_grid2 <- prior_variance * c(0.1, 0.5, 1.0, 1.5, 2.0)
  uvbsers <- purrr::map(sigma2_grid2, ~fit_uvb_ser(sim$X, sim$y, prior_variance=.x, estimate_prior_variance=F))
  log_bfs <- do.call('cbind', purrr::map(uvbsers, ~.x$lbf))
  log_bf_ser <- do.call('cbind', purrr::map(uvbsers, ~.x$lbf_model))
  # here we plot the in the laplace approximait
  par(mfrow=c(1,2))
  plot(sigma2_grid2, log_bfs[1,], ylim = c(-1, 1) + range(log_bfs), type = 'l',  ylab = 'UVB log(BF)', xlab = 'prior variance', main='Individual variables')
  for(i in 2:nrow(log_bfs)){
    lines(sigma2_grid2, log_bfs[i,], ylim = c(-1, 1) + range(log_bfs), type = 'l')
  }
  lines(sigma2_grid2, log_bf_ser, col='red')
  plot(sigma2_grid2, log_bf_ser, col='red', type='l', ylab = 'UVB log(BF_SER)', xlab = 'prior variance', main='SER')


  # QUADRATURE
  # lets look at how Gauss-Hermite quadrature performs at a range of values (0, 1, 10, 66.39, 120)
  x <- sim$X[,66]
  y <- sim$y
  jj_fit <- fit_univariate_vb(x, y, tau0 = 1/prior_variance)
  elbo <- tail(jj_fit$elbos, 1)
  b0 <- jj_fit$delta
  mu <- jj_fit$mu
  var <- 1/jj_fit$tau

  #' do Gauss-Hermite quadrature for a range of # of quadrature points
  #' centered on the posterior mean, s is a rescaling of the posterior variance
  quad_scale_variance <- function(ns, s, mu, var){
    purrr::map_dbl(ns, ~compute_log_marginal(
      x, y, b0 = b0, mu = mu, var = s*var, prior_variance = prior_variance, n=.x, verbose = F))
  }

  # suggests we can improve quadrature performance by
  s1 <- quad_scale_variance(ns, 1, mu, var)
  s2 <- quad_scale_variance(ns, 2, mu, var)
  s4 <- quad_scale_variance(ns, 4, mu, var)
  s8 <- quad_scale_variance(ns, 8, mu, var)
  s16 <- quad_scale_variance(ns, 16, mu, var)

  ll0 <- sum(dbinom(y, 1, mean(y), log=T))
  ylim <- range(c(elbo - ll0, s1 - ll0, s2 - ll0, s4 - ll0, s8 - ll0, s16 - ll0))


  plot(log2(ns), s1 - ll0, type='l', col='red', ylim=ylim, ylab = 'log p(y)')
  lines(log2(ns), s2 - ll0, col='blue')
  lines(log2(ns), s4 - ll0, col='green')
  lines(log2(ns), s8 - ll0, col='purple')
  lines(log2(ns), s16 - ll0, col='orange')
  abline(h = elbo - ll0, lty=2) # plot elbo
  abline(h = f3(prior_variance))
  abline(h = f5(prior_variance))

  par(mfrow = c(1, 1))
  quadser <- with(sim, fit_quad_ser(X, y, prior_variance=prior_variance, glm_ser=uvbser, n=16))
  plot(quadser$lbf, uvbser$lbf, main='Quad vs Univarite JJ bound'); abline(0, 1, col='red')

  abline(h = elbo, lty=2) # plot elbo

  legend(x ="bottomright",
         legend=c("s=1", "s=2", "s=4", "s=8", "s=16", 'elbo'),
         col = c("red","blue", "green", "purple", "orange","black"),
         lty = c(1, 1, 1, 1, 1, 2))

  # a
  log_joint <- function(b){
    psi <- b0 + x * b
    return(sum(y * psi - log(1 + exp(psi))) + dnorm(b, 0, sqrt(prior_variance), log=T))
  }
  # a

  xi <- jj_fit$xi
  vb_bound <- function(b){
    psi <- b0 + x * b
    xs <- 2 * (y - 0.5) * psi
    lambda_xi <- (0.5 - sigmoid(psi)) / (2 * xi)
    bound <- sum(log(sigmoid(xi)) + 0.5 * (xs - xi) + lambda_xi *(xs^2 - xi^2)) + dnorm(b, 0, sqrt(prior_variance), log=T)
    return(bound)
  }

  bs <- seq(mu-1, mu+0.5, by=0.01)

  mu_opt <- optimize(log_joint, c(-10, 10), maximum = T)$maximum
  plot(bs, Vectorize(log_joint)(bs), type='l')
  abline(v = mu_opt, lty=3)

  lines(bs, Vectorize(vb_bound)(bs), type='l', col='blue')
  abline(v = mu, lty=3, col='blue')
  legend(
    x = 'topright',
    legend = c('log p(y, b)', 'JJ approximation', 'Posterior mode', 'JJ posterior mean'),
    col = c('black', 'blue'),
    lty = c(1, 1, 3, 3),
    cex=0.5
  )

  # suggests we can improve quadrature performance by
  par(mfrow=c(1,1))
  s1 <- quad_scale_variance(ns, 1, mu_opt, var)
  s2 <- quad_scale_variance(ns, 2, mu_opt, var)
  s4 <- quad_scale_variance(ns, 4, mu_opt, var)
  s8 <- quad_scale_variance(ns, 8, mu_opt, var)
  s16 <- quad_scale_variance(ns, 16, mu_opt, var)

  ylim <- range(c(elbo, s1, s2, s4, s8, s16))

  plot(log2(ns), s1, type='l', col='red', ylim=ylim, ylab = 'log p(y)')
  lines(log2(ns), s2, col='blue')
  lines(log2(ns), s4, col='green')
  lines(log2(ns), s8, col='purple')
  lines(log2(ns), s16, col='orange')
  abline(h = elbo, lty=2) # plot elbo

  legend(x ="bottomright",
         legend=c("s=1", "s=2", "s=4", "s=8", "s=16", 'elbo'),
         col = c("red","blue", "green", "purple", "orange","black"),
         lty = c(1, 1, 1, 1, 1, 2))

  optimize(log_joint, c(-10, 10), maximum = T)

  f2 <- function(var0){with(glmser, compute_log_labf_ser(betahat, shat2, lr, var0, 1/length(betahat)))}
  with(sim, fit_uvb_ser(X, y, prior_variance=prior_variance, estimate_prior_variance=T))
})

test_that("Search for residual variance UVB-SER", {
  library(tictoc)

  sim <- sim_ser()

  # estimate prior variance with UVB
  f <- function(prior_variance){
    with(sim, fit_uvb_ser(X, y, prior_variance=prior_variance, estimate_prior_variance=F))$lbf_model
  }
  prior_variance <- optimize(f, c(1e-10, 100), maximum = T)$maximum
  plot(prior_variance*seq(0.5, 2, 0.25), Vectorize(f)(prior_variance*seq(0.5, 2, 0.25)))

  tic()
  uvbser <- with(sim, fit_uvb_ser(X, y, prior_variance=prior_variance, estimate_prior_variance=F))
  toc()

  # estimate prior variance with Laplace MLE
  glmser <- with(sim, fit_glm_ser(X, y, prior_variance=1, estimate_prior_variance=T))
  abfser <- with(sim, fit_glm_ser(X, y, prior_variance=1, estimate_prior_variance=T, laplace = F))

  glue::glue('UVB estimated prior variance = {uvbser$prior_variance}\n Laplace estimated prior variance = {glmser$prior_variance}\n Wakefield estimated prior variance = {abfser$prior_variance}')

  # LAPLACE
  f2 <- function(var0){with(glmser, compute_log_labf(betahat, shat2, lr, var0))}
  sigma2_grid <- seq(0.1, 50, by=0.1)
  laplace_log_bfs <- do.call('cbind', purrr::map(sigma2_grid, f2))

  f3 <- function(var0){with(glmser, compute_log_labf_ser(betahat, shat2, lr, var0, pi = 1/length(betahat)))}
  laplace_log_bf_sers <- purrr::map_dbl(sigma2_grid, f3)

  # here we plot the in the laplace approximait
  par(mfrow=c(1,2))
  plot(sigma2_grid, laplace_log_bfs[1,], ylim = c(-1, 1) + range(laplace_log_bfs), type = 'l',  ylab = 'Laplace log(BF)', xlab = 'prior variance', main='Individual variables')
  for(i in 2:nrow(laplace_log_bfs)){
    lines(sigma2_grid, laplace_log_bfs[i,], ylim = c(-1, 1) + range(laplace_log_bfs), type = 'l')
  }
  lines(sigma2_grid, laplace_log_bf_sers, col='red')
  plot(sigma2_grid, laplace_log_bf_sers, col='red', type='l', ylab = 'Laplace log(BF_SER)', xlab = 'prior variance', main='SER')


  # ABF
  f2 <- function(var0){with(glmser, compute_log_abf(betahat, shat2, var0))}
  sigma2_grid <- seq(0.1, 50, by=0.1)
  log_abfs <- do.call('cbind', purrr::map(sigma2_grid, f2))

  f3 <- function(var0){with(glmser, compute_log_abf_ser(betahat, shat2, var0, pi = 1/length(betahat)))}
  log_abf_ser <- purrr::map_dbl(sigma2_grid, f3)

  # here we plot the in the laplace approximait
  par(mfrow=c(1,2))
  plot(sigma2_grid, laplace_log_bfs[1,], ylim = c(-1, 1) + range(log_abfs), type = 'l',  ylab = 'Wakefield log(BF)', xlab = 'prior variance', main='Individual variables')
  for(i in 2:nrow(log_abfs)){
    lines(sigma2_grid, log_abfs[i,], ylim = c(-1, 1) + range(log_abfs), type = 'l')
  }
  lines(sigma2_grid, log_abf_ser, col='red')
  plot(sigma2_grid, log_abf_ser, col='red', type='l', ylab = 'Wakefield log(BF_SER)', xlab = 'prior variance', main='SER')


  # QUADRATURE
  # lets look at how Gauss-Hermite quadrature performs at a range of values (0, 1, 10, 66.39, 120)
  x <- sim$X[,1]
  y <- sim$y
  jj_fit <- fit_univariate_vb(x, y, 0, 0, 0, 1, tau0 = 1/prior_variance)
  elbo <- tail(jj_fit$elbos, 1)
  b0 <- jj_fit$delta
  mu <- jj_fit$mu
  var <- 1/jj_fit$tau

  #' do Gauss-Hermite quadrature for a range of # of quadrature points
  #' centered on the posterior mean, s is a rescaling of the posterior variance
  quad_scale_variance <- function(ns, s, mu, var){
    purrr::map_dbl(ns, ~compute_log_marginal(
      x, y, b0 = b0, mu = mu, var = s*var, prior_variance = prior_variance, n=.x, verbose = F))
  }

  # suggests we can improve quadrature performance by
  s1 <- quad_scale_variance(ns, 1, mu, var)
  s2 <- quad_scale_variance(ns, 2, mu, var)
  s4 <- quad_scale_variance(ns, 4, mu, var)
  s8 <- quad_scale_variance(ns, 8, mu, var)
  s16 <- quad_scale_variance(ns, 16, mu, var)

  ylim <- range(c(elbo, s1, s2, s4, s8, s16))

  plot(log2(ns), s1, type='l', col='red', ylim=ylim, ylab = 'log p(y)')
  lines(log2(ns), s2, col='blue')
  lines(log2(ns), s4, col='green')
  lines(log2(ns), s8, col='purple')
  lines(log2(ns), s16, col='orange')
  abline(h = elbo, lty=2) # plot elbo

  legend(x ="bottomright",
         legend=c("s=1", "s=2", "s=4", "s=8", "s=16", 'elbo'),
         col = c("red","blue", "green", "purple", "orange","black"),
         lty = c(1, 1, 1, 1, 1, 2))


  # a
  log_joint <- function(b){
    psi <- b0 + x * b
    return(sum(y * psi - log(1 + exp(psi))) + dnorm(b, 0, sqrt(prior_variance), log=T))
  }
  bs <- seq(mu-1, mu+0.5, by=0.1)
  plot(bs, Vectorize(log_joint)(bs), type='l')
  mu_opt <- optimize(log_joint, c(-10, 10), maximum = T)$maximum
  abline(v = mu)
  abline(v = mu_opt)

  # suggests we can improve quadrature performance by
  par(mfrow=c(1,1))
  s1 <- quad_scale_variance(ns, 1, mu_opt, var)
  s2 <- quad_scale_variance(ns, 2, mu_opt, var)
  s4 <- quad_scale_variance(ns, 4, mu_opt, var)
  s8 <- quad_scale_variance(ns, 8, mu_opt, var)
  s16 <- quad_scale_variance(ns, 16, mu_opt, var)

  ylim <- range(c(elbo, s1, s2, s4, s8, s16))

  plot(log2(ns), s1, type='l', col='red', ylim=ylim, ylab = 'log p(y)')
  lines(log2(ns), s2, col='blue')
  lines(log2(ns), s4, col='green')
  lines(log2(ns), s8, col='purple')
  lines(log2(ns), s16, col='orange')
  abline(h = elbo, lty=2) # plot elbo

  legend(x ="bottomright",
         legend=c("s=1", "s=2", "s=4", "s=8", "s=16", 'elbo'),
         col = c("red","blue", "green", "purple", "orange","black"),
         lty = c(1, 1, 1, 1, 1, 2))

  optimize(log_joint, c(-10, 10), maximum = T)

  f2 <- function(var0){with(glmser, compute_log_labf_ser(betahat, shat2, lr, var0, 1/length(betahat)))}
  with(sim, fit_uvb_ser(X, y, prior_variance=prior_variance, estimate_prior_variance=T))

  sim <- make_c1_sim()
  glmser <- with(sim, logisticsusie::fit_glm_ser(X, y, prior_variance=1))
  uvbser <- with(sim, logisticsusie::fit_uvb_ser(X, y, prior_variance=1))
  quadser <- with(sim, logisticsusie::fit_quad_ser(X, y, prior_variance = 1, glm_ser=uvbser, n=16))
})

