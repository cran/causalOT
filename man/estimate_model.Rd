% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treatment_effect.R
\name{estimate_model}
\alias{estimate_model}
\title{Function to estimate outcome models}
\usage{
estimate_model(data, causalWeights, model.function, separate.estimation, ...)
}
\arguments{
\item{data}{A \code{\link[=dataHolder]{dataHolder()}} object}

\item{causalWeights}{A \linkS4class{causalWeights} object}

\item{model.function}{The model function passed by the user}

\item{separate.estimation}{TRUE or FALSE, should models be estimated separately in each group?}

\item{...}{Extra agruments passed to the predict functions}
}
\value{
a list with slots \code{y_hat_0}, \code{y_hat_1}, and \code{fit}.
}
\description{
Function to estimate outcome models
}
\keyword{internal}
