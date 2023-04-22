#' Sigmoid
#' sigmoid function coverts log-odds to probability
#' @param x Log odds
#' @return Returns the probability
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' @export
logodds <- function(p) {
  return(log(p) - log(1 - p))
}

logSumExp <- matrixStats::logSumExp

softmax <- function(x) {
  return(exp(x - logSumExp(x)))
}

.monotone <- function(v) {
  return(all(tail(v, -1) - head(v, -1) >= 0))
}

rowCumSum <- function(X) {
  do.call(rbind, apply(X, 1, cumsum, simplify = F))
}

.clamp <- function(v, .upper = 100, .lower = -100) {
  return(pmax(.lower, pmin(v, .upper)))
}

is.odd <- function(x) {
  x %% 2 == 1
}
is.even <- function(x) {
  x %% 2 == 0
}

