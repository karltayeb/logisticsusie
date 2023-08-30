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
  diff <- Inf
  # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
  while (diff > tol) {

    beta_post_old <- beta_post

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

    # diff to monitor convergence
    diff <- mean((beta_post - beta_post_old)**2)

    if (iter > maxit) {
      warning("Maximum number of iterations reached")
      break
    }
  }
  timer <- tictoc::toc()

  # now, get intercept w/ MLE, holding our final estimate of beta to be fixed
  beta <- colSums(alpha * mu)
  pred <- (X %*% beta)[, 1]

  # TODO: get intercept for whole mode

  # add the final ser fits
  names(fits) <- paste0("L", 1:L)
  prior_variance <- purrr::map_dbl(fits, ~purrr::pluck(.x, 'prior_variance'))

  # compute PIP using only effects with prior variance != 0
  if(sum(prior_variance > 0) == 0){
    pip <- rep(0, p)
  } else{
    pip <- compute_pip(alpha[prior_variance > 0,])
  }

  # compute cs
  cs <- compute_cs(alpha)

  res <- list(
    alpha = alpha,
    mu = mu,
    var = var,
    prior_variance = prior_variance,
    fits = fits,
    iter = iter,
    elapsed_time = unname(timer$toc - timer$tic),
    q_history = head(q_history, iter),
    pip = pip,
    cs = cs
  )
  return(res)
}

#' Generalized IBSS
#'
#' approximate GLM SuSiE using generalized IBSS heuristic
#' by default this is for logistic regression but you can specify
#' what GLM to use via "family" argument
#' dots get passed to `ibss_from_ser` see documentation for options
#' @param X design matrix
#' @param y response
#' @param L number of single effects
#' @param laplace boolean to use Laplace approximation to BF rather than ABF--
#'  we recommend keeping set to default `TRUE`
#' @param estimate_prior_variance boolean to estimate prior variance
#' @param family family for glm
#' @export
generalized_ibss <- function(X, y, L=10, laplace=T, estimate_prior_variance=T, family='binomial', ...){
  ser_fun <- purrr::partial(fit_glm_ser2,
                            laplace=laplace,
                            estimate_prior_variance=estimate_prior_variance,
                            family=family)
  ibss_from_ser(X, y, L=L, ser_function = ser_fun, ...)
}
