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

get_cs <- function(alpha, requested_coverage = 0.95) {
  rho <- order(-alpha)
  idx <- min(sum(cumsum(alpha[rho]) < requested_coverage) + 1, length(alpha))
  cs <- rho[1:idx]
  coverage <- sum(alpha[cs])
  return(list(cs = cs, prob = alpha[rho[1:idx]], size = idx, requested_coverage = requested_coverage, coverage = coverage))
}

get_all_cs <- function(alpha, requested_coverage = 0.95) {
  L <- dim(alpha)[1]
  sets <- purrr::map(1:L, ~ get_cs(alpha[.x, ], requested_coverage))
  names(sets) <- paste0("L", 1:L)
  return(sets)
}

# convenient table of CSs from get_all_cs2
cs_tbl2 <- function(alpha) {
  get_all_cs(alpha) %>%
    dplyr::tibble() %>%
    tidyr::unnest_wider(1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(top_feaure = cs[1], top_feature_alpha = prob[1]) # %>% unnest_longer(c(cs, prob))
}


