#' Initialize q with lasso
#'
#' Providies initialization for q using lasso
initialize_q_with_lasso <- function(X, y, L='auto', family='binomial'){
  if(is.null(colnames(X))){
    colnames(X) <- paste0('var', 1:ncol(X))
  }

  lasso_fit <- glmnet::cv.glmnet(X, y, family=family)
  Coef <- coef(lasso_fit, s=lasso_fit$lambda) # matrix of coefficients along regularization path
  coef_path <- rev(sort(rowSums(Coef[2:nrow(Coef),] != 0))) # order that variables enter the model
  if(L == 'auto'){
    # use the "1se" penalty
    lasso_coef <- tail(coef(lasso_fit, s='lambda.1se'), -1)[,1]
    include <- names(which(lasso_coef[names(coef_path)] != 0))
    mu_lasso <- lasso_coef[include]
    L <- length(include)
  } else {
    n_coef <- purrr::map_int(1:ncol(Coef), ~sum(Coef[,.x]!=0)) - 1
    lambda <- lasso_fit$lambda[which(cumsum(n_coef >= L) == 1)]
    lasso_coef <- tail(coef(lasso_fit, s=lambda)[,1], -1)

    # add effects in order that they enter lasso
    include <- names(which(lasso_coef[names(coef_path)] != 0)[1:L])
    mu_lasso <- lasso_coef[include]
  }

  idx <- match(include, colnames(X))
  p <- ncol(X)
  mu <- matrix(data = 0, nrow = L, ncol = p)
  mu[1:length(idx), idx] <- diag(mu_lasso)

  # variance doesn't matter for G-IBSS, except for estimating lbf
  var <- matrix(data = 1, nrow = L, ncol = p)

  alpha <- matrix(data = 1/p, nrow=L, ncol=p)
  alpha[1:length(idx),] <- 0
  alpha[1:length(idx), idx] <- diag(rep(1, length(idx)))

  q <- list(mu=mu, var=var, alpha=alpha, lasso_fit = lasso_fit)
  return(q)
}

#' @export
gibss_lasso_init <- function(X, y, L = 10, family = 'binomial', ...){
  message('initializing with lasso...')
  q_init <- initialize_q_with_lasso(X, y, L, family)
  L <- nrow(q_init$alpha) # in case L = 'auto', infer from output of q_init

  fit <- gibss(X, y, L, family=family, init = q_init, ...)
  return(fit)
}

#' @export
gibss <- function(X, y, L = 10,
                  prior_variance=1,
                  estimate_prior_variance=T,
                  family = 'binomial',
                  augment = T,
                  laplace = T,
                  glm_mapper = map_fastgml,
                  maxit=100,
                  init = NULL,
                  tol=1e-8){
  message('fitting GLM-SuSiE via generalized IBSS...')
  ser_fun <- purrr::partial(fit_glm_ser2,
                            estimate_prior_variance = estimate_prior_variance,
                            family = family,
                            laplace = laplace,
                            augment = augment)
  fit <- ibss_from_ser(X, y, L, ser_function = ser_fun, init = init, maxit=maxit, tol=tol)
  fit$lasso_fit <- q_init$lasso_fit  # return lasso fit too
  fit$cs <- compute_cs(fit$alpha)
  fit$prior_variance <- purrr::map_dbl(fit$fits, 'prior_variance')
  fit$lbf_ser <- purrr::map_dbl(fit$fits, 'lbf_model')
  fit$cs_size <- purrr::map_int(fit$cs, ~length(.x$cs))
  return(fit)
}
