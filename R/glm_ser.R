# Adapted from: https://andrewg3311.github.io/susieR_logistic_wflow/susie_logistic_demonstration.html#adjustments_to_susie_code

#' Map univariate regression
#'
#' Fit univariate glm for each column of X
#' @param X n x p design matrix
#' @param y n vector response
#' @param o optional offset, vector of length n
#' @param estimate_intercept boolean to fit intercept
#' @param family family for glm, e.g. binomial
#' @returns list containing intercept and effect estimates,
#'    standard errors, and likelihood ratios
map_univariate_regression <- function(X, y, off = NULL, estimate_intercept=T, family='binomial', augment=F){
  if(augment==T){
    warning('`augment == T` not implemeneted')
  }
  p <- ncol(X)
  betahat <- numeric(p)
  shat2 <- numeric(p)
  intercept <- rep(0, p)
  lr <- rep(0, p)

  if (is.null(off)) {
    off <- rep(0, length(y))
  }

  for (j in 1:p) {
    if (estimate_intercept) {
      log.fit <- glm(y ~ X[, j] + 1 + offset(off), family = family)
      intercept[j] <- unname(coef(log.fit)[1])
    }
    else {
      log.fit <- glm(y ~ X[, j] - 1 + offset(off), family = family)
    }
    log.fit.coef <- summary(log.fit)$coefficients
    betahat[j] <- ifelse(is.na(log.fit.coef[1 + estimate_intercept,
                                            1]), 0, log.fit.coef[1 + estimate_intercept, 1])
    shat2[j] <- ifelse(is.na(log.fit.coef[1 + estimate_intercept,
                                          2]), Inf, log.fit.coef[1 + estimate_intercept, 2]^2)
    lr[j] <- 0.5 * (log.fit$null.deviance - log.fit$deviance)
  }

  res <- list(
    p = p,
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept,
    lr = lr
  )
  return(res)
}


#' Prepare data for fast glm
#'
#' Helper function for preparing the design matrix, response vector, and offset
#' @param x the variable of interest
#' @param y the response vector
#' @param o fixed offset, length(o) == length(y)
#' @param intercept boolean to include column of 1s in design matrix
#' @param augment boolean to augment data to ensure existence of MLE
#' @returns a list containing design `X`, reponse `y`, and offset `o`
prep_data_for_fastglm <- function(x, y, o, intercept=T, augment=F){
  if(augment){
    n <- length(y) + 4
    x <- c(x, c(1, 1, 0, 0))
    y <- c(y, c(1, 0, 1, 0))
    o <- c(o, rep(0, 4))
    weights <- c(rep(1, n-4), rep(1e-6, 4))
  } else{
    n <- length(y)
    weights <- rep(1, n)
  }

  if(intercept){
    ones <- rep(1, n)
    X <- cbind(ones, x)
  } else{
    X <- as.matrix(x, nrow=n)
  }
  regression_data <- list(X=X, y=y, o=o, weights=weights)
  return(regression_data)
}

#' Fit fast GLM
#'
#' Fit fast GLM, and extract necessary information for SER
#' @param x the variable of interest
#' @param y the response vector
#' @param o fixed offset, length(o) == length(y)
#' @param family the family of the glm, e.g. "binomial", "poisson"
#' @param augment boolean to augment data to ensure existence of MLE
#' @param fit0
fit_fast_glm <- function(x, y, o, ll0, family='binomial'){
  # 1. fit model
  fit <- fastglm::fastglm(y = y,
                          x = as.matrix(cbind(rep(1, length(y)), x)),
                          offset = o,
                          family = family)
  # 2. processing
  coef <- unname(fit$coefficients)
  intercept <- coef[1]
  betahat <- coef[2]
  shat2 <-  fit$se[2]^2
  ll <- -0.5 * fit$deviance
  lr  <- ll - ll0
  res <- list(
    betahat = betahat,
    shat2 = shat2,
    intercept = intercept,
    lr = lr,
    ll0 = ll0,
    ll = ll
  )
  return(res)
}

#' Map univariate regression fit with fastglm
#'
#' Fit univariate glm for each column of X
#' @param X n x p design matrix
#' @param y n vector response
#' @param o optional offset, vector of length n
#' @param estimate_intercept boolean to fit intercept (not implemented)
#' @param family family for glm, e.g. binomial
#' @returns list containing intercept and effect estimates,
#'    standard errors, and likelihood ratios
map_fastglm <- function(X, y, o = rep(0, length(y)), estimate_intercept=T, family='binomial'){
  if(!estimate_intercept){
    warning("estimate_intercept = FALSE not implemented")
  }
  p <- ncol(X)

  fit0 <- fastglm::fastglm(y = y,
                          x = matrix(rep(1, length(y)), ncol=1),
                          offset = o,
                          family = family)
  ll0 <- -0.5 * fit0$deviance

  # 1. fit each variable
  res <- purrr::map_dfr(1:p, ~fit_fast_glm(X[, .x], y, o, ll0, family=family))
  res$p <- p
  return(res)
}

#' compute SER using GLM, use Laplace approximation to BF rather than ABF
#' @export
fit_glm_ser <- function(X, y, o = NULL,
                         prior_variance = 1, intercept = T,
                         prior_weights = NULL, family = "binomial",
                         laplace = T,
                         estimate_prior_variance=T,
                         min_prior_variance = 1e-3,
                         glm_mapper = map_fastglm) {

  # univariate regressions
  uni <- glm_mapper(X, y, o,
                    estimate_intercept=intercept,
                    family=family)
  lr <- uni$lr
  betahat <- uni$betahat
  shat2 <- uni$shat2
  intercept <- uni$intercept
  p <- ncol(X)

  # default pi
  if (is.null(prior_weights)) {
    prior_weights <- rep(1/p, p)
  }

  res <- asymptotic_ser(betahat, shat2, intercept, lr, prior_weights, estimate_prior_variance=estimate_prior_variance, laplace=T, min_prior_variance=min_prior_variance)
  class(res) <- 'asymptotic_ser'
  return(res)
}
