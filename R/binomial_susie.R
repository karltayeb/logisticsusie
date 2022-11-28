
###
# Computations ----
##


#' Get Binomial SER coefficients
#' Extract regression coefficients from binomial SER fit
#' @param fit Binomial SER object
#' @return Return E[\beta]
coef.binsusie <- function(fit, k, idx = NULL) {
  if (missing(k)) {
    b <- colSums(.get_alpha(fit, idx) * .get_mu(fit, idx))
  } else {
    b <- colSums(.get_alpha(fit$logreg_list[[k]], idx) * .get_mu(fit$logreg_list[[k]], idx))
  }

  return(b)
}

#' Expected linear prediction
compute_Xb.binsusie <- function(fit, k) {
  if (missing(k)) {
    Xb <- drop(fit$data$X %*% colSums(.get_alpha(fit) * .get_mu(fit)))
    Zd <- drop(fit$data$Z %*% colSums(.get_delta(fit)))
  } else {
    Xb <- drop(fit$data$X %*% colSums(.get_alpha(fit$logreg_list[[k]]) * .get_mu(fit$logreg_list[[k]])))
    Zd <- drop(fit$data$Z %*% colSums(.get_delta(fit$logreg_list[[k]])))
  }

  return(Matrix::drop(Xb + Zd))
}


compute_Xb2.binsusie <- function(fit, k) {
  if (missing(k)) {
    B <- .get_alpha(fit) * .get_mu(fit)
    XB <- fit$data$X %*% t(B)
    Xb <- rowSums(XB)
    Zd <- drop(fit$data$Z %*% colSums(.get_delta(fit)))

    B2 <- .get_alpha(fit) * (.get_mu(fit)^2 + .get_var(fit))
    b2 <- colSums(B2)

    Xb2 <- fit$data$X2 %*% b2 + Xb**2 - rowSums(XB^2)
    Xb2 <- drop(Xb2 + 2 * Xb * Zd + Zd^2)
  } else {
    B <- .get_alpha(fit$logreg_list[[k]]) * .get_mu(fit$logreg_list[[k]])
    XB <- fit$data$X %*% t(B)
    Xb <- rowSums(XB)
    Zd <- drop(fit$data$Z %*% colSums(.get_delta(fit$logreg_list[[k]])))

    B2 <- .get_alpha(fit$logreg_list[[k]]) * (.get_mu(fit$logreg_list[[k]])^2 + .get_var(fit$logreg_list[[k]]))
    b2 <- colSums(B2)

    Xb2 <- fit$data$X2 %*% b2 + Xb**2 - rowSums(XB^2)
    Xb2 <- drop(Xb2 + 2 * Xb * Zd + Zd^2)
  }

  return(Xb2)
}


#' Compute KL(q(\beta) || p(\beta)) for Sum of SERs
#' Note this does not include KL(q(\omega) || p(\omega))
#' since that gets computes as part of the JJ bound
compute_kl.binsusie <- function(fit, k) {
  if (missing(k)) {
    kl <- sum(purrr::map_dbl(seq(fit$hypers$L), ~ .compute_ser_kl(
      .get_alpha(fit, .x),
      .get_pi(fit, .x),
      .get_mu(fit, .x),
      .get_var(fit, .x),
      .get_prior_mean(fit, .x),
      .get_var0(fit, .x)
    )))
  } else {
    kl <- sum(purrr::map_dbl(seq(fit$logreg_list[[k]]$hypers$L), ~ .compute_ser_kl(
      .get_alpha(fit$logreg_list[[k]], .x),
      .get_pi(fit$logreg_list[[k]], .x),
      .get_mu(fit$logreg_list[[k]], .x),
      .get_var(fit$logreg_list[[k]], .x),
      .get_prior_mean(fit$logreg_list[[k]], .x),
      .get_var0(fit$logreg_list[[k]], .x)
    )))
  }

  return(kl)
}


#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
#' NOTE: this only works when xi is updated!
jj_bound.binsusie <- function(fit, k) {
  if (missing(k)) {
    xi <- .get_xi(fit)
    Xb <- compute_Xb.binsusie(fit)
    kappa <- compute_kappa(fit)
    n <- .get_N(fit)

    bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi)
  } else {
    xi <- .get_xi(fit$logreg_list[[k]])
    Xb <- compute_Xb.binsusie(fit, k = k)
    kappa <- compute_kappa(fit, k = k)
    n <- .get_N(fit)

    bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi)
  }
  return(bound)
}

compute_elbo.binsusie <- function(fit, k) {
  if (missing(k)) {
    jj <- sum(jj_bound.binsusie(fit)) # E[log p(y, w | b) - logp(w)]
    kl <- compute_kl.binsusie(fit) # KL(q(b) || p(b))
  } else {
    jj <- sum(jj_bound.binsusie(fit, k = k)) # E[log p(y, w | b) - logp(w)]
    kl <- compute_kl.binsusie(fit, k = k) # KL(q(b) || p(b))
  }

  return(jj - kl)
}

#' Compute E[p(z, w| b) - q(w)], which is the same as
#' the bound on the logistic function proposed by Jaakkola and Jordan
jj_bound2.binsusie <- function(fit, k) {
  if (missing(k)) {
    xi <- .get_xi(fit)
    Xb <- compute_Xb.binsusie(fit)
    kappa <- compute_kappa(fit)
    n <- .get_N(fit)

    Xb2 <- compute_Xb2.binsusie(fit)
    omega <- compute_omega(fit)

    bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi) + 0.5 * omega * (xi^2 - Xb2)
  } else {
    xi <- .get_xi(fit$logreg_list[[k]])
    Xb <- compute_Xb.binsusie(fit, k = k)
    kappa <- compute_kappa(fit, k = k)
    n <- .get_N(fit)

    Xb2 <- compute_Xb2.binsusie(fit, k = k)
    omega <- compute_omega(fit, k = k)

    bound <- n * log(sigmoid(xi)) + (kappa * Xb) - (0.5 * n * xi) + 0.5 * omega * (xi^2 - Xb2)
  }


  return(bound)
}

compute_elbo2.binsusie <- function(fit, k) {
  if (missing(k)) {
    jj <- sum(jj_bound2.binsusie(fit)) # E[log p(y, w | b) - logp(w)]
    kl <- compute_kl.binsusie(fit) # KL(q(b) || p(b))
  } else {
    jj <- sum(jj_bound2.binsusie(fit, k = k)) # E[log p(y, w | b) - logp(w)]
    kl <- compute_kl.binsusie(fit, k = k) # KL(q(b) || p(b))
  }

  return(jj - kl)
}

###
# Updates ----
###


#' update alpha, mu, and var
update_b.binsusie <- function(fit, k, fit_intercept = TRUE, fit_prior_variance = TRUE, fit_alpha = TRUE, update_idx = NULL) {
  if (missing(k)) {
    if (is.null(update_idx)) {
      update_idx <- seq(fit$hypers$L)
    }
    shift <- compute_Xb.binsusie(fit)
    for (l in update_idx) {
      # remove current effect estimate
      shift <- Matrix::drop(shift - compute_Xb.binser(fit, idx = l))

      # update SER
      post_l <- update_b.binser(fit, idx = l, shift = shift)
      fit$params$mu[l, ] <- post_l$mu
      fit$params$var[l, ] <- post_l$var

      if (fit_alpha) {
        fit$params$alpha[l, ] <- post_l$alpha
      }

      # update intercept/fixed effect covariates
      if (fit_intercept) {
        fit$params$delta[l, ] <- update_delta.binser(fit, idx = l, shift = shift)
      }

      # update prior_variance
      if (fit_prior_variance) {
        fit$hypers$prior_variance[l] <- update_prior_variance.binser(fit, idx = l)
      }

      # add current effect estimate
      shift <- shift + compute_Xb.binser(fit, k = k, idx = l)
    }
  } else {
    if (is.null(update_idx)) {
      update_idx <- seq(fit$logreg_list[[k]]$hypers$L)
    }
    shift <- compute_Xb.binsusie(fit, k = k)
    for (l in update_idx) {
      # remove current effect estimate
      shift <- Matrix::drop(shift - compute_Xb.binser(fit, k = k, idx = l))

      # update SER
      post_l <- update_b.binser(fit, k = k, idx = l, shift = shift)
      fit$params$mu[l, ] <- post_l$mu
      fit$params$var[l, ] <- post_l$var

      if (fit_alpha) {
        fit$params$alpha[l, ] <- post_l$alpha
      }

      # update intercept/fixed effect covariates
      if (fit_intercept) {
        fit$params$delta[l, ] <- update_delta.binser(fit, k = k, idx = l, shift = shift)
      }

      # update prior_variance
      if (fit_prior_variance) {
        fit$hypers$prior_variance[l] <- update_prior_variance.binser(fit, k = k, idx = l)
      }

      # add current effect estimate
      shift <- shift + compute_Xb.binser(fit, k = k, idx = l)
    }
  }

  return(fit)
}


#' update for variational parameter parameter xi q(w) = PG(N, xi)
update_xi.binsusie <- function(fit, k) {
  if (missing(k)) {
    Xb2 <- compute_Xb2.binsusie(fit)
    xi <- sqrt(abs(Xb2))
  } else {
    Xb2 <- compute_Xb2.binsusie(fit, k = k)
    xi <- sqrt(abs(Xb2))
  }

  return(xi)
}

###
# Fitting ----
###

.init.binsusie.params <- function(n, p, p2, L) {
  params <- list(
    alpha = matrix(rep(1, p * L) / p, nrow = L),
    mu = matrix(rep(0, L * p), nrow = L),
    var = matrix(rep(1, L * p), nrow = L),
    delta = matrix(rep(0, p2 * L, nrow = L)),
    xi = rep(1e-3, n),
    tau = 0 # initialized at end of this function
  )
  return(params)
}

.init.binsusie.hypers <- function(n, p, L, prior_mean, prior_variance, prior_weights) {
  stopifnot("`prior_mean` must be scalar of length L`" = length(prior_mean) %in% c(1, L))
  if (length(prior_mean) == 1) {
    prior_mean <- rep(prior_mean, L)
  }
  stopifnot("`prior_variance` must be scalar of length L`" = length(prior_variance) %in% c(1, L))
  if (length(prior_variance) == 1) {
    prior_variance <- rep(prior_variance, L)
  }

  if (is.null(prior_weights)) {
    prior_weights <- matrix(rep(1 / p, p * L), nrow = L)
  } else if (is.matrix(prior_weights)) {
    stopifnot("`prior weights must be a L x p matrix or p vector" = all(dim(prior_weights) == c(L, p)))
    stopifnot("`prior weights must sum to one" = all(rowSums(prior_weights) == 1))
  } else {
    stopifnot("`prior weights must be a L x p matrix or p vector" = length(prior_weights) == p)
    stopifnot("`prior weights must sum to one" = sum(prior_weights) == 1)
    prior_weights <- matrix(rep(prior_weights, L), nrow = L, byrow = T)
  }

  hypers <- list(
    L = L,
    pi = prior_weights, # categorical probability of nonzero effect
    prior_mean = prior_mean, # prior mean of effect
    prior_variance = prior_variance # prior variance of effect, TODO: include as init arg,
  )
  return(hypers)
}

#' initialize SER
init.binsusie <- function(data, L = 5, prior_mean = 0, prior_variance = 1, prior_weights = NULL, kidx = NULL) {
  n <- nrow(data$X)
  p <- ncol(data$X)
  p2 <- ncol(data$Z)

  params <- .init.binsusie.params(n, p, p2, L)
  hypers <- .init.binsusie.hypers(n, p, L, prior_mean, prior_variance, prior_weights)

  data$X2 <- data$X^2

  # TODO: check that data has y, N, X, Z
  fit.init <- list(
    data = data,
    params = params,
    hypers = hypers,
    elbo = c(-Inf)
  )
  fit.init$kidx <- kidx
  fit.init$params$tau <- compute_tau(fit.init)
  return(fit.init)
}


iter.binsusie <- function(fit, k, fit_intercept = TRUE, fit_prior_variance = TRUE, fit_xi = TRUE, fit_alpha = TRUE) {
  if (missing(k)) {
    # update b
    fit <- update_b.binsusie(fit,
      fit_intercept = fit_intercept,
      fit_prior_variance = fit_prior_variance,
      fit_alpha = fit_alpha
    )

    # update xi
    if (fit_xi) {
      fit$params$xi <- update_xi.binsusie(fit)
      fit$params$tau <- compute_tau(fit)
    }
  } else {
    # update b
    fit$logreg_list[[k]] <- update_b.binsusie(fit,
      k                  = k,
      fit_intercept      = fit_intercept,
      fit_prior_variance = fit_prior_variance,
      fit_alpha          = fit_alpha
    )

    # update xi
    if (fit_xi) {
      fit$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit, k)
      fit$logreg_list[[k]]$params$tau <- compute_tau(fit, k)
    }
  }


  return(fit)
}

#' Fit the binomial single effect regression
fit.binsusie <- function(fit,
                         k,
                         maxiter = 10,
                         tol = 1e-3,
                         fit_intercept = TRUE,
                         fit_prior_variance = TRUE,
                         fit_xi = TRUE,
                         fit_alpha = TRUE,
                         fast_elbo = TRUE,
                         kidx = NULL) {
  if (missing(k)) {
    for (i in 1:maxiter) {
      fit <- iter.binsusie(fit,
        fit_intercept      = fit_intercept,
        fit_prior_variance = fit_prior_variance,
        fit_xi             = fit_xi,
        fit_alpha          = fit_alpha
      )
      # update elbo
      if (fast_elbo) {
        fit$elbo <- c(fit$elbo, compute_elbo.binsusie(fit))
      } else {
        fit$elbo <- c(fit$elbo, compute_elbo2.binsusie(fit))
      }
      if (.converged(fit, tol)) {
        break
      }
    }
  } else {
    for (i in 1:maxiter) {
      fit$logreg_list[[k]] <- iter.binsusie(fit,
        k                  = k,
        fit_intercept      = fit_intercept,
        fit_prior_variance = fit_prior_variance,
        fit_xi             = fit_xi,
        fit_alpha          = fit_alpha
      )
      # update elbo
      if (fast_elbo) {
        fit$logreg_list[[k]]$elbo <- c(fit$logreg_list[[k]]$elbo, compute_elbo.binsusie(fit, k))
      } else {
        fit$logreg_list[[k]]$elbo <- c(fit$logreg_list[[k]]$elbo, compute_elbo2.binsusie(fit, k))
      }
      if (.converged(fit, tol)) {
        break
      }
    }
  }

  return(fit)
}

#' remove a single effect and refit model
#' This function zeros out a single effect and refits the model
#' Used as a subroutine for `prune_model` which removes irrelevant features
prune_model_idx <- function(fit, k, idx, fit_intercept = T, fit_prior_variance = F, tol = 1e-3, maxiter = 10) {
  if (missing(k)) {
    # copy the susie model
    fit2 <- fit

    # don't update this index
    update_idx <- seq(fit2$hypers$L)
    update_idx <- update_idx[idx != update_idx]

    # instead, fix in to 0
    fit2$params$mu[idx, ] <- 0
    fit2$params$var[idx, ] <- 1e-10
    fit2$params$delta[idx, ] <- 0
    fit2$hypers$prior_mean[idx] <- 0
    fit2$hypers$prior_variance[idx] <- 1e-10

    # new elbo
    fit2$elbo <- compute_elbo.binsusie(fit2)

    # CAVI loop
    for (i in 1:maxiter) {
      # update b
      fit2 <- update_b.binsusie(fit2,
        fit_intercept = fit_intercept,
        fit_prior_variance = fit_prior_variance,
        update_idx = update_idx
      )
      # update xi
      fit2$params$xi <- update_xi.binsusie(fit2)
      fit2$params$tau <- compute_tau(fit2)

      # update elbo
      fit2$elbo <- c(fit2$elbo, compute_elbo.binsusie(fit2))
      if (.converged(fit2, tol)) {
        break
      }
    }
  } else {
    # copy the susie model
    fit2 <- fit
    fit$logreg_list
    # don't update this index
    update_idx <- seq(fit2$hypers$L)
    update_idx <- update_idx[idx != update_idx]

    # instead, fix in to 0
    fit2$logreg_list[[k]]$params$mu[idx, ] <- 0
    fit2$logreg_list[[k]]$params$var[idx, ] <- 1e-10
    fit2$logreg_list[[k]]$params$delta[idx, ] <- 0
    fit2$logreg_list[[k]]$hypers$prior_mean[idx] <- 0
    fit2$logreg_list[[k]]$hypers$prior_variance[idx] <- 1e-10

    # new elbo
    fit2$logreg_list[[k]]$elbo <- compute_elbo.binsusie(fit2)

    # CAVI loop
    for (i in 1:maxiter) {
      # update b
      fit2$logreg_list[[k]] <- update_b.binsusie(fit2,
        k = k,
        fit_intercept = fit_intercept,
        fit_prior_variance = fit_prior_variance,
        update_idx = update_idx
      )
      # update xi
      fit2$logreg_list[[k]]$params$xi <- update_xi.binsusie(fit2, k)
      fit2$logreg_list[[k]]$params$tau <- compute_tau(fit2, k)

      # update elbo
      fit2$logreg_list[[k]]$elbo <- c(fit2$logreg_list[[k]]$elbo, compute_elbo.binsusie(fit2, k))
      if (.converged(fit2$logreg_list[[k]], tol)) {
        break
      }
    }
  }

  return(fit2)
}

#' Prune Model
#' Remove irrelevant features from logistic SuSiE fit
#' @param fit a fit logistic Susie model
#' @param check_null_threshold if ELBO_null + check_null_threshold > ELBO drop the single effect
#' @param fit_intercept boolean to re-estimate intercept in simplified model, defaul TRUE
#' @param fit_prior_variance boolean to re-estimate prior variance in simplified model, default FALSE
#' @param tol tolerance to declare convergence when fitting simplified model
prune_model <- function(fit, k, check_null_threshold, fit_intercept = T, fit_prior_variance = F, tol = 1e-3) {
  if (missing(k)) {
    # order by explained variance
    var0 <- purrr::map_dbl(seq(fit$hypers$L), ~ update_prior_variance.binser(fit, idx = .x))
    idx_order <- order(var0, decreasing = F)

    for (idx in idx_order) {
      fit2 <- prune_model_idx(fit, idx = idx, fit_intercept = fit_intercept, fit_prior_variance = fit_prior_variance, tol = tol)
      null_diff <- tail(fit2$elbo, 1) - tail(fit$elbo, 1)
      if (null_diff + check_null_threshold > 0) {
        message(paste0("Null ELBO diff: ", null_diff, "... REMOVING effect"))
        fit <- fit2
      } else {
        message(paste0("Null ELBO diff: ", null_diff, "... KEEPING effect"))
      }
    }
  } else {
    # order by explained variance
    var0 <- purrr::map_dbl(
      seq(fit$logreg_list[[k]]$hypers$L),
      ~ update_prior_variance.binser(fit, k = k, idx = .x)
    )
    idx_order <- order(var0, decreasing = F)

    for (idx in idx_order) {
      fit2 <- prune_model_idx(fit,
        k = k,
        idx = idx,
        fit_intercept = fit_intercept,
        fit_prior_variance = fit_prior_variance,
        tol = tol
      )
      null_diff <- tail(fit2$elbo, 1) - tail(fit$elbo, 1)
      if (null_diff + check_null_threshold > 0) {
        message(paste0("Null ELBO diff: ", null_diff, "... REMOVING effect"))
        fit <- fit2
      } else {
        message(paste0("Null ELBO diff: ", null_diff, "... KEEPING effect"))
      }
    }
  }

  return(fit)
}


#' BinSuSiE Wrapup
#' Compute fit summaries (PIPs, CSs, etc) and
#' Organize binsusie fit so that it is compatable with `susieR` function
binsusie_wrapup <- function(fit, k, prior_tol = 0) {
  if (missing(k)) {
    class(fit) <- c("binsusie", "susie")

    # TODO put back into original X scale
    X <- fit$data$X

    fit$alpha <- fit$params$alpha
    fit$mu <- fit$params$mu
    fit$mu2 <- fit$params$mu^2 + fit$params$var
    fit$V <- fit$hypers$prior_variance

    colnames(fit$alpha) <- colnames(X)
    colnames(fit$mu) <- colnames(X)

    fit$pip <- susieR::susie_get_pip(fit, prior_tol = prior_tol)
    names(fit$pip) <- colnames(X)

    fit$sets <- susieR::susie_get_cs(fit, X = X)
    fit$intercept <- colSums(fit$params$delta)[1]
  } else {
    class(fit) <- c("binsusie", "susie")

    # TODO put back into original X scale
    X <- fit$data$X

    fit$alpha <- fit$logreg_list[[k]]$params$alpha
    fit$mu <- fit$logreg_list[[k]]$params$mu
    fit$mu2 <- fit$logreg_list[[k]]$params$mu^2 + fit$logreg_list[[k]]$params$var
    fit$V <- fit$logreg_list[[k]]$hypers$prior_variance

    colnames(fit$logreg_list[[k]]$alpha) <- colnames(X)
    colnames(fit$logreg_list[[k]]$mu) <- colnames(X)

    fit$logreg_list[[k]]$pip <- susieR::susie_get_pip(fit$logreg_list[[k]],
      prior_tol = prior_tol
    )
    names(fit$logreg_list[[k]]$pip) <- colnames(X)

    fit$logreg_list[[k]]$sets <- susieR::susie_get_cs(fit$logreg_list[[k]], X = X)
    fit$logreg_list[[k]]$intercept <- colSums(fit$logreg_list[[k]]$params$delta)[1]
  }




  return(fit)
}

# help -----

#' @export
binsusie_get_pip <- function(fit, k, prune_by_cs = FALSE, prior_tol = 1e-09) {
  if (missing(k)) {
    # TODO: filter out components that don't need to be included (see Gao's suggestion)
    include_idx <- which(fit$hypers$prior_variance > prior_tol)
    if (length(include_idx) > 0) {
      pip <- susieR::susie_get_pip(fit$params$alpha[include_idx, , drop = F])
    } else {
      pip <- rep(0, dim(fit$params$alpha)[2])
    }
  } else {
    # TODO: filter out components that don't need to be included (see Gao's suggestion)
    include_idx <- which(fit$logreg_list[[k]]$hypers$prior_variance > prior_tol)
    if (length(include_idx) > 0) {
      pip <- susieR::susie_get_pip(fit$logreg_list[[k]]$params$alpha[include_idx, , drop = F])
    } else {
      pip <- rep(0, dim(fit$logreg_list[[k]]$params$alpha)[2])
    }
  }

  return(pip)
}

#' @export
binsusie_plot <- function(fit, k, y = "PIP") {
  if (missing(k)) {
    res <- with(fit, list(alpha = params$alpha, pip = pip, sets = sets))
    class(res) <- "susie"
    susieR::susie_plot(res, y)
  } else {
    res <- with(
      fit$logreg_list[[k]],
      list(alpha = params$alpha, pip = pip, sets = sets)
    )
    class(res) <- "susie"
    susieR::susie_plot(res, y)
  }
}

#' Binomial SuSiE
#' Fit Binomial SuSiE via coordinate ascent variational inference
#' @param X a n x p matrix of covariates
#' @param y an n vector of integer counts, bernoulli/binomial observations
#' @param N the number of binomial trials, defaults to 1, may be a scalar or vector of length n
#' @param Z fixed effect covaraites (including intercept). If null just a n x 1 matrix of ones
#' @param scale if TRUE, scale the columns of X to unit variate
#' @param center if TRUE, center the columns of X to mean zero
#' @param prior_mean the prior mean of each non-zero element of b. Either a scalar or vector of length L.
#' @param prior_variance the prior variance of each non-zero element of b. Either a scalar or vector of length L. If `estimate_prior_variance=TRUE` the value provides an initial estimate of prior variances
#' @param prior_weights prior probability of selecting each column of X, vector of length p summing to one, or an L x p matrix
#' @param intercept
#' @param estimate_prior_variance
#' @param s_init a logistic susie object to initialize with, NOTE if non-null, we ignore `prior_mean`, `prior_variance`, and `prior_weights`
#' @param returns a fit Binomial SuSiE model, which is compatable with summary functions from `susieR` package
#' @export
binsusie <- function(X,
                     y,
                     N = rep(1, length(y)), # number of trials for each
                     Z = NULL,
                     L = min(10, ncol(X)),
                     scale = FALSE,
                     center = TRUE,
                     prior_mean = 0.0, # set prior mean (feature added for Gao and Bohan)
                     prior_variance = 1.0, # TODO: make scaled prior variance?
                     prior_weights = NULL, # vector of length `p` gjvj g the prior probability of each column of X having nonzero effect... = hypers$pi
                     # null_weight = NULL,
                     # standardize = TRUE,
                     intercept = TRUE,
                     estimate_prior_variance = TRUE,
                     # estimate_prior_method = c("optim", "EM", "simple"),  right now only EM is implimented
                     check_null_threshold = 0,
                     prior_tol = 1e-09,
                     prune = FALSE,
                     s_init = NULL, # previous fit with which to initialize NO
                     coverage = 0.95,
                     min_abs_corr = 0.5,
                     max_iter = 100,
                     tol = 0.001,
                     verbose = FALSE,
                     # track_fit = FALSE,
                     # refine = FALSE,
                     n_purity = 100) {
  # Initialize model
  if (is.null(s_init)) {
    fit <- binsusie_init(X, y, N, Z, L, scale, center, prior_mean, prior_variance, prior_weights)
  } else {
    fit <- s_init
  }

  # model fitting
  for (i in 1:max_iter) {
    fit <- iter.binsusie(
      fit,
      fit_intercept = intercept,
      fit_prior_variance = estimate_prior_variance
    )
    fit$elbo <- c(fit$elbo, compute_elbo.binsusie(fit))
    if (.converged(fit, tol)) {
      break
    }
  }

  # model pruning-- removes irrelevant components
  if (prune) {
    fit <- prune_model(fit, check_null_threshold, intercept, estimate_prior_variance, tol = tol)
  }

  # wrap up (computing PIPs, CSs, etc)
  fit <- binsusie_wrapup(fit, prior_tol)
  return(fit)
}
