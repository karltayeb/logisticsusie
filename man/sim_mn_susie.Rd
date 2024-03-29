% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{sim_mn_susie}
\alias{sim_mn_susie}
\title{Simulate from a multinomial logit model.
logits for each category are sparse linear functions of X (L non-zero effects per category)
multinomial sample taken with probability of softmax(logits)}
\usage{
sim_mn_susie(n = 1000, p = 50, L = 3, N = 1, K = 10)
}
\description{
Simulate from a multinomial logit model.
logits for each category are sparse linear functions of X (L non-zero effects per category)
multinomial sample taken with probability of softmax(logits)
}
