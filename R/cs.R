compute_cs1 <- function (alpha, requested_coverage = 0.95) {
  rho <- order(-alpha)
  idx <- min(sum(cumsum(alpha[rho]) < requested_coverage) +
               1, length(alpha))
  cs <- rho[1:idx]
  coverage <- sum(alpha[cs])
  res <- list(cs = cs, alpha = alpha[rho[1:idx]], size = idx,
              requested_coverage = requested_coverage, coverage = coverage)
  class(res) <- 'cs'
  return(res)
}

#' @export
print.cs <- function(cs){
  cs_str <- paste(head(cs$cs, 6), collapse=', ')
  if(length(cs$cs) > 6){
    cs_str <- paste0(cs_str, ', ...')
  }
  ln1 <- paste0('CS = {', cs_str, '}')
  alpha_str <- paste(head(round(cs$alpha, 2), 6), collapse=', ')
  if(length(cs$cs) > 6){
    alpha_str <- paste0(alpha_str, ', ...')
  }
  ln2 <- paste0('alpha = {', alpha_str, '}')
  ln3 <- paste0(
    'size = ', cs$size,
    ', coverage = ', round(cs$coverage, 3),
    ', requested = ', round(cs$requested_coverage, 3))
  cat(ln1, '\n', ln2, '\n', ln3)
}

#' Compute credible sets
#'
#' Compute credible sets at a requested coverage level
#' @param alpha a vector or matrix, where each row is a probability vector
#' @param requested_coverage the minimum probability mass the CS must satisfy
#' @returns a credible set object with the credible set and other summaries
#'     if `alpha` is a vector, returns a list of credible set objects with names
#'     `L1`, `L2`, ... etc.
#' @export
compute_cs <- function (alpha, requested_coverage = 0.95) {
  if(is.null(dim(alpha))){
    return(compute_cs1(alpha, requested_coverage = requested_coverage))
  } else{
    L <- dim(alpha)[1]
    sets <- purrr::map(1:L, ~compute_cs(alpha[.x, ], requested_coverage))
    names(sets) <- paste0("L", 1:L)
    return(sets)
  }
}

# convenient table of CSs from get_all_cs2
cs_tbl2 <- function(alpha) {
  compute_cs(alpha) %>%
    dplyr::tibble() %>%
    tidyr::unnest_wider(1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(top_feaure = cs[1], top_feature_alpha = prob[1]) # %>% unnest_longer(c(cs, prob))
}

