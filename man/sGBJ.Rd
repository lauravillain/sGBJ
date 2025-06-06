% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sGBJ.R
\name{sGBJ}
\alias{sGBJ}
\title{Compute the sGBJ statistic and its p-value quantifying a gene set expression association
with survival}
\usage{
sGBJ(surv, factor_matrix, covariates = NULL, nperm = 300)
}
\arguments{
\item{surv}{a \code{\link[survival]{Surv}} object of length \code{n}}

\item{factor_matrix}{a \code{n x p} \code{data.frame} of the expression for the
particular gene set of interest being tested}

\item{covariates}{a \code{n x l} matrix of the covariates to adjust upon. Default is \code{NULL}}

\item{nperm}{number of permutations performed to estimate the \code{epsilon} matrix.
Default is \code{300}.}
}
\value{
The sGBJ statistic and its associated p-value associated
}
\description{
This function is the main function of the sGBJ package to perform Gene Set Analysis
in the context of time-to-event outcome.
}
\examples{
n <- 100
surv_data <- data.frame(Time = runif(n = n, min = 0, max = 100),
                        event = rbinom(n = n, size = 1, prob = 0.5))
surv <- survival::Surv(time = surv_data$Time, event = surv_data$event)

factor_matrix <- data.frame(P1 = rnorm(n = n),
                             P2 = rnorm(n = n))

sGBJ::sGBJ(surv,factor_matrix, nperm = 2)

}
