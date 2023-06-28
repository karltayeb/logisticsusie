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
ibss_from_ser <- function(X, y, L = 10, prior_variance = 1., prior_weights = NULL, tol = 1e-3, maxit = 100, estimate_intercept = TRUE, ser_function = NULL, init = NULL) {
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
  beta_post_history <- vector(mode = "list", length = maxit)
  iter <- 0
  # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
  while (ibss_l2(beta_post_history, iter - 1) > tol) {
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
    beta_post_history[[iter]] <- beta_post

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
    beta_post_history = head(beta_post_history, iter)
  )
  return(res)
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
