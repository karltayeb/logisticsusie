prepare_data <- function(x, ...){
  UseMethod("prepare_data", x)
}

.check_X <- function(X) {
  stopifnot("data$X must be matrix/Matrix" = inherits(X, c("matrix", "Matrix")))
}

#' @export
binsusie_prep_data <- function(X, y, N, Z, shift=0, shift_var=0, center = TRUE, scale = FALSE) {
  # center and scale data
  .check_X(X) # throws errors if something is wrong
  X <- Matrix::Matrix(scale(X, center = center, scale = scale))

  # If N was passed as a scalar, convert to vector of length n
  n <- length(y)
  if (length(N) == 1) {
    N <- rep(N, n)
  }

  # set Z to intercept covariate if not provided
  if (is.null(Z)) {
    Z <- matrix(rep(1, n), nrow = n)
  }

  # Make data list
  # TODO: store X means and standard errors
  data <- list(X = X, X2 = X^2, Z = Z, y = y, N = N, shift=shift, shift_var=shift_var)
  return(data)
}


#' Number of binomial trials for a stick breaking representation of multinomial
#' @param Y an N x K matrix of assignment probabilities
#' @return an N x K matrix
stick_breaking_N <- function(Y) {
  cumY <- do.call(rbind, apply(Y, 1, cumsum, simplify = F))
  K <- ncol(Y)
  Nk <- cumY[, K] - cumY + Y
  return(Nk)
}

#' Prepare data for multinomial regression
#' @param X matrix of covoariates
#' @param y matrix of multinomial observations (each row an observation)
#' @param N vector giving number of trials for each multinomial observation
#' @param Z matrix of fixed-effect covariates (including intercept)
#' @params shift, shift_var mean and varianance respectively of individual-specific random intercept
#' @param center boolean to center `X`
#' @param scale boolean to scale `X` to unit variance
mnsusie_prep_data <- function(X, Y, Z=NULL, shift=0, shift_var=0, center = TRUE, scale = FALSE) {
  # center and scale data
  .check_X(X) # throws errors if something is wrong
  X <- Matrix::Matrix(scale(X, center = center, scale = scale))

  # If N was passed as a scalar, convert to vector of length n
  n <- nrow(Y)

  # set Z to intercept covariate if not provided
  if (is.null(Z)) {
    Z <- matrix(rep(1, n), nrow = n)
  }

  Nk <- stick_breaking_N(Y)

  # Make data list
  # TODO: store X means and standard errors
  data <- list(X = X, X2 = X^2, Z = Z, Y = Y, Nk = Nk, shift=0, shift_var=0)
  return(data)
}
