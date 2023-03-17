prepare_data <- function(x, ...){
  UseMethod("prepare_data", x)
}

.check_X <- function(X) {
  stopifnot("data$X must be matrix/Matrix" = inherits(X, c("matrix", "Matrix")))
}

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


