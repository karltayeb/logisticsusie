% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ser_from_univariate.R
\name{compute_log_abf_ser}
\alias{compute_log_abf_ser}
\title{Compute log of the asymptotic Bayes factor (ABF) for the SER vs null model}
\usage{
compute_log_abf_ser(betahat, s2, prior_variance, pi)
}
\arguments{
\item{betahat}{vector of MLEs for the effect size}

\item{s2}{vector of standard errors}

\item{prior_variance}{variance of the effect size b ~ N(0, prior_variance)}

\item{pi}{prior over non-zero effects}
}
\value{
log of the ABF for the SER log(\sum_j \pi_j ABF_j)
}
\description{
Compute log of the asymptotic Bayes factor (ABF) for the SER vs null model
}
