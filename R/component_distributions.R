# Generics -----

#' Convolved logpdf
#'
#' Generic function for `component_distribution` objects, computes the loglikelihood of an observation from
#' the component distribution, corrupted by gaussian measurement error
#'
#' @param betahat observations
#' @param se standard errors of observations
#' @return the convolved log density log p(betahat | se) = log \int N(betahat | beta, se) p(beta) dbeta
convolved_logpdf <- function(x, ...) {
  UseMethod("convolved_logpdf", x)
}

#' Update Params
#'
#' Generic function for updating the parameters of a `component_distribution`
#'
#' @param betahat observations
#' @param se standard errors of observations
#' @param weights weights to weigh each observation by when optimizing
#' @return a new `component_distribution` object with update parameters
update_params <- function(x, ...) {
  UseMethod("update_params", x)
}

# Defaults ----

convolved_logpdf.component_distribution <- function(dist, betahat, se) {
  # TODO: compute convolved pdf by numerical integration?
  # Will probably require logpdf/pdf be implimented for each component
  stop("Generic convolved logpdf not implimented yet")
}

update_params.component_distribution <- function(dist, betahat, se) {
  # TODO: compute convolved pdf by numerical integration?
  # Will probably require logpdf/pdf be implimented for each component
  stop("Generic update not implimented")
}

# Point mass Component----

point_component <- function(mu = 0) {
  f <- list(mu = mu)
  class(f) <- c("point", "component_distribution")
  return(f)
}

is.point <- function(x) {
  inherits(x, "point")
}

convolved_logpdf.point <- function(dist, betahat, se) {
  # return(dnorm(betahat, sd=se, log=T))
  sd <- se
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 1e3, -1e3)
  return(logp)
}

# Normal Component----

normal_component <- function(mu = 0, var = 1) {
  f <- list(mu = mu, var = var)
  class(f) <- c("normal", "component_distribution")
  return(f)
}

is.normal <- function(x) {
  inherits(x, "normal")
}

convolved_logpdf.normal <- function(dist, betahat, se) {
  sd <- sqrt(se^2 + dist$var)
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 100, -100)
  return(logp)
}

update_params.normal <- function(dist, betahat, se, weights) {
  # TODO: right now it's just a grid search, but we can impliment EM update easily
  var.init <- dist$var
  grid_var <- 2^(seq(-1, 1, by = 0.1) + log2(var.init))
  grid_dist <- purrr::map(grid_var, ~ normal_component(mu = dist$mu, var = .x))
  ll <- purrr::map_dbl(grid_dist, ~ sum(convolved_logpdf.normal(.x, betahat, se) * weights, na.rm = T))
  return(grid_dist[[which.max(ll)]])
}



# beta components
beta_component <- function(alpha = 0.5) {
  f <- list(alpha = alpha)
  class(f) <- c("beta", "component_distribution")
  return(f)
}


convolved_logpdf.beta <- function(dist, p, se = 1) {
  logp <- dbeta(p, 1, dist$alpha)
  logp <- .clamp(logp, 100, -100)
  return(logp)
}

update_params.beta <- function(dist, p, se = 1, weights) {
  # TODO: decide how to update alpha
  # weights come from posterior assignment probabilities.
}
