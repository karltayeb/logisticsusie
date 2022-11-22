
# Not in operator

'%!in%' <- function(x,y)!('%in%'(x,y))

#' Sigmoid
#' sigmoid function coverts log-odds to probability
#' @param x Log odds
#' @return Returns the probability
sigmoid <- function(x) {
  return(1 / (1 + exp(-x)))
}

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

#' helper function gets difference in ELBO between iterations
.diff <- function(fit) {
  return(diff(tail(fit$elbo, 2)))
}

#' check if ELBO is converged to some tolerance threshold
.converged <- function(fit, tol = 1e-3) {
  is.converged <- F
  if ((length(fit$elbo) > 1) & .diff(fit) <= tol) {
    is.converged <- T
  }
  return(is.converged)
}


rowCumSum <- function(X) {
  do.call(rbind, apply(X, 1, cumsum, simplify = F))
}


#' Get Credible Sets
#' Wraps susieR::susie_get_cs
#' @export
binsusie_get_cs <- function(fit,
                            coverage = 0.95,
                            min_abs_corr = 0.5,
                            dedup = TRUE,
                            squared = FALSE,
                            check_symmetric = TRUE,
                            n_purity = 100) {
  res <- list(
    X = fit$data$X,
    alpha = fit$params$alpha,
    V = fit$hypers$prior_variance
  )
  cs <- susieR::susie_get_cs(res,
    X = fit$data$X,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    dedup = dedup,
    squared = squared,
    check_symmetric = check_symmetric,
    n_purity = n_purity
  )
  return(cs)
}

# implement normal and point mass component distributions

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

get_all_cs <- function(fit, requested_coverage = 0.95) {
  sets <- purrr::map(1:fit$hypers$L, ~ get_cs(fit$params$alpha[.x, ], requested_coverage))
  names(sets) <- paste0("L", 1:fit$hypers$L)
  return(sets)
}







#' @title Preparing data for mococomo fit
#' @details Preparing data for mococomo fit for two type of input p-values or estimated regression coefficients
#' currently only support one type of entry, whether betahat or p-values
#'
#' @param betahat numeric, vector of regression coefficients  list of containing the following element see
#' @param se numeric, corresponding standard error
#' @param p numeric observed p-values
#' @param maxiter numeric, maximum numerous of iteration set to 100 by defaults
#' @param tol tolerance in term of change in ELBO value for stopping criterion
#' @export
#' @example
#' see \link{\code{fit.mococomo}}
set_data_mococomo <- function (betahat, se,p, X, ... ){

  if( !is.numeric( X))
  {
    stop("X should be numercial vector")
  }
  if( !is.matrix(X))
  {
    stop("X should be a matrix")
  }


  if( !missing(betahat)){
    if( !is.numeric( betahat))
    {
      stop("betahat should be numercial vector")
    }
    if( !is.numeric( se))
    {
      stop("se should be numercial vector")
    }
    if( ! (sum(c(length(se)==length(betahat), length(se)==nrow(X)))==2  ))
    {
      stop(" The number of lines in X should be equal to the number of entries in Betahat and se ")
    }
    dat <-  list(betahat= betahat,
                 se=se ,
                 X=X)
    class(dat) <- c("normal","data_mococomo")
  }


  if(!missing(p)){
    if( !is.numeric( p))
    {
      stop("p should be numercial vector")
    }
    if( !  (length(p)==nrow(X))==1   )
    {
      stop(" The number of lines in X should be equal to the number of entries in p ")
    }
    dat <-  list(p= p,
                 se=rep(1, length(p)) ,
                 X=X)
    class(dat) <- c("beta","data_mococomo")
  }


  return(dat)
}

