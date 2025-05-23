% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper-functions.R
\name{compute_weights}
\alias{compute_weights}
\title{Compute Weights for Newton Step}
\usage{
compute_weights(p, epsilon_small = 1e-08)
}
\arguments{
\item{p}{Predicted probabilities.}

\item{epsilon_small}{A threshold to avoid weights that are too small.}
}
\value{
A vector of weights.
}
\description{
This function computes the weights used in the Newton-Raphson update
for logistic regression.
}
\keyword{internal}
