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
ibss_from_ser <- function(X, y, L = 10, tol = 1e-8, maxit = 100, ser_function = NULL, init = NULL) {
  if (is.null(ser_function)) {
    stop("You need to specify a fit function `fit_glm_ser`, `fit_vb_ser`, etc")
  }

  p <- ncol(X)
  n <- nrow(X)

  if(is.null(init)){
    init <- null_initialize_ibss(n, p, L, 1.)
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
    if (iter > maxit) {
      warning("Maximum number of iterations reached")
      break
    }

    beta_post_old <- beta_post

    for (l in 1:L) {
      # remove effect from previous iteration
      fixed <- fixed - (X %*% beta_post[l, ])[, 1]

      # fit SER
      ser_l <- ser_function(X, y, o = fixed)

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
    diff <- max((beta_post - beta_post_old)**2)
  }
  timer <- tictoc::toc()

  # now, get intercept w/ MLE, holding our final estimate of beta to be fixed
  beta <- colSums(alpha * mu)
  pred <- (X %*% beta)[, 1]

  # TODO: properly estimate intercept of whole model
  intercept <- with(ser_l, sum(intercept * alpha))

  # add the final ser fits
  names(fits) <- paste0("L", 1:L)
  prior_variance <- purrr::map_dbl(fits, ~purrr::pluck(.x, 'prior_variance'))
  lbf_ser <- purrr::map_dbl(fits, ~purrr::pluck(.x, 'lbf_model'))

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
    intercept = intercept,
    prior_variance = prior_variance,
    fits = fits,
    iter = iter,
    elapsed_time = unname(timer$toc - timer$tic),
    q_history = head(q_history, iter),
    pip = pip,
    lbf_ser = lbf_ser,
    cs = cs
  )
  class(res) <- 'generalized_ibss'
  return(res)
}


#' @export
predict.generalized_ibss <- function(fit, X){
  psi <- with(fit, intercept + drop(X %*% colSums(alpha * mu)))[, 1]
  return(psi)
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
generalized_ibss <- function(X, y, L=10, tol=1e-8, maxit=100, init=NULL, ...){
  # make SER function for GLM, uses asymptotic approximation
  ser_fun <- purrr::partial(fit_glm_ser, ...)

  # fit IBSS using the SER function
  ibss_from_ser(X, y, L=L, ser_function = ser_fun, tol = tol, maxit = maxit, init=init)
}


#' Generalized IBSS Quadrature
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
generalized_ibss_quad <- function(X, y, L=10, tol=1e-8, maxit=100, init=NULL, ...){
  # make SER function for GLM, uses asymptotic approximation
  ser_fun <- purrr::partial(fit_quad_ser, ...)

  # fit IBSS using the SER function
  ibss_from_ser(X, y, L=L, ser_function = ser_fun, tol = tol, maxit = maxit, init=init)
}
