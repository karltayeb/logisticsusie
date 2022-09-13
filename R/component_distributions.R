# Generics -----
convolved_logpdf <- function(x, ...) {
  UseMethod("convolved_logpdf", x)
}


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

#' Convolved log density
#' @param  dist is a list with parameters for the normal component
#' @param betahat vector of oberveations
#' @param se vector of standard errors of observations
#' @return A vector of log densities log p(\hat\beta | se) = log N \hat\beta ; \mu, sigma^2 + se^2)
convolved_logpdf.normal <- function(dist, betahat, se) {
  sd <- sqrt(se^2 + dist$var)
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 100, -100)
  return(logp)
}


#' Update variance of the normal parameter on a grid centered around the current estimate
#' on a log2 scale
update_params.normal <- function(dist, betahat, se, weights) {
  var.init <- dist$var
  grid_var <- 2^(seq(-1, 1, by = 0.1) + log2(var.init))
  grid_dist <- purrr::map(grid_var, ~ normal_component(mu = dist$mu, var = .x))
  ll <- purrr::map_dbl(grid_dist, ~ sum(convolved_logpdf.normal(.x, betahat, se) * weights, na.rm = T))
  return(grid_dist[[which.max(ll)]])
}
