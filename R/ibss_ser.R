# Functions for a generic IBSS algorithm
# allows us to specify an SER function to run IBSS with...

null_initialize_ibss <- function(n, p, L, prior_variance){
  # place to store posterior info for each l = 1, ..., L
  alpha <- matrix(1/p, nrow = L, ncol = p)
  mu <- matrix(0, nrow = L, ncol = p)
  var <- matrix(prior_variance, nrow = L, ncol = p)
  init <- list(alpha = alpha, mu = mu, var = var)
}

#' IBSS from SER
#'
#' Fit IBSS for aribtrary SER, using expected predictions as a fixed offset
#' @param X design matrix
#' @param y response
#' @param L number of effects
#' @param prior_variance variance of each single effect, if `estimate_prior_variance = T` this is an initialization
#' @param prior_weights the prior probability of selecting each variable (columns of `X`)
#' @param tol tolerance to declare convergence
#' @param maxit maximum number of iterations
#' @param estimate_intercept boolean to estimate intercept
#' @param ser_fun function for performing SER, must return PIPs, `alpha` and posterior mean `mu`
#' @param init initialization
#' @export
ibss_from_ser <- function(X, y, L = 10, prior_variance = 1., prior_weights = NULL, tol = 1e-8, maxit = 100, estimate_intercept = TRUE, ser_function = NULL, init = NULL) {
  if (is.null(ser_function)) {
    stop("You need to specify a fit function `fit_glm_ser`, `fit_vb_ser`, etc")
  }

  p <- ncol(X)
  n <- nrow(X)

  if(is.null(init)){
    init <- null_initialize_ibss(n, p, L, prior_variance)
  }

  alpha <- init$alpha
  mu <- init$mu
  var <- init$var

  beta_post <- mu * alpha  # expected posterior effect
  fixed <- (X %*% colSums(beta_post))[, 1] # expected predictions from SERs

  tictoc::tic() # start timer
  fits <- vector(mode = "list", length = L)
  q_history <- vector(mode = "list", length = maxit + 1)
  q_history[[1]] <- list(alpha = alpha, mu = mu, var = var)

  iter <- 1
  kl_diff <- Inf
  # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
  while (kl_diff > tol) {
    for (l in 1:L) {
      # remove effect from previous iteration
      fixed <- fixed - (X %*% beta_post[l, ])[, 1]

      # fit SER
      ser_l <- ser_function(X, y,
        o = fixed,
        prior_variance = prior_variance,
        intercept = estimate_intercept,
        prior_weights = prior_weights
      )

      # store
      alpha[l, ] <- ser_l$alpha
      mu[l, ] <- ser_l$mu
      var[l, ] <- ser_l$var

      # update beta_post
      beta_post[l, ] <- ser_l$alpha * ser_l$mu

      # add back new fixed portion
      fixed <- fixed + (X %*% beta_post[l, ])[, 1]

      # save current ser
      fits[[l]] <- ser_l
    }

    iter <- iter + 1
    q_history[[iter]] <- list(alpha = alpha, mu = mu, var = var)

    # compute kl between q after this iter, and q before this iter
    kl_diff <- compute_sum_ser_kls(q_history[[iter]], q_history[[iter - 1]])

    if (iter > maxit) {
      warning("Maximum number of iterations reached")
      break
    }
  }
  timer <- tictoc::toc()

  # now, get intercept w/ MLE, holding our final estimate of beta to be fixed
  beta <- colSums(alpha * mu)
  pred <- (X %*% beta)[, 1]
  int <- coef(glm(y ~ 1 + offset(pred), family = "binomial"))

  # add the final ser fits
  names(fits) <- paste0("L", 1:L)

  res <- list(
    alpha = alpha,
    mu = mu,
    var = var,
    intercept = int,
    fits = fits,
    iter = iter,
    elapsed_time = unname(timer$toc - timer$tic),
    q_history = head(q_history, iter)
  )
  return(res)
}


compute_sum_ser_kls <- function(q1, q2){
  L <- nrow(q1$alpha)
  kl <- 0
  for(l in 1:L){
    kl12 <- compute_ser_kl(
      q1$alpha[l,], q2$alpha[l,], q1$mu[l,], q1$var[l,], q2$mu[l,], q2$var[l,])
    kl21 <- compute_ser_kl(
      q2$alpha[l,], q1$alpha[l,], q2$mu[l,], q2$var[l,], q1$mu[l,], q1$var[l,])
    kl <- kl + 0.5 * (kl12 + kl21)
  }
  return(kl)
}

ibss_l2 <- function(beta_post_history, iter) {
  if (iter < 2) {
    res <- Inf
  } else {
    res <- norm(
      beta_post_history[[iter - 1]] - beta_post_history[[iter]],
      type = "2"
    )
  }
  return(res)
}

#' @export
ibss_monitor_convergence <- function(fit) {
  purrr::map_dbl(2:(fit$iter), ~ ibss_l2(fit$beta_post_history, .x))
}
