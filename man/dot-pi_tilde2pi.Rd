% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_function.R
\name{.pi_tilde2pi}
\alias{.pi_tilde2pi}
\title{Convert \tilde pi to \pi
Converts conditional probabilities p(z = k | z >=k) k = 1... K
to marginal probilities p(z =k) k = 1... K}
\usage{
.pi_tilde2pi(pi_tilde)
}
\arguments{
\item{pi_tilde}{vector of conditional probabilities (except for last one p(z=K | Z >=K) = 1)}
}
\value{
pi a vector of probabilities summing to one
}
\description{
Convert \tilde pi to \pi
Converts conditional probabilities p(z = k | z >=k) k = 1... K
to marginal probilities p(z =k) k = 1... K
}
