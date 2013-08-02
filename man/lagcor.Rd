\name{lagcor}
\alias{lagcor}
\title{Lagged pearson}
\usage{lagcor(x, y, lags=-2:2, ...)}
\description{(best) Lagged pearson}
\value{the best correlation over the considered lags}
\author{henning}
\arguments{\item{x}{a numeric vector}
\item{y}{a numeric vector}
\item{lags}{lags to consider for the x-vector}
\item{...}{passed on to cor}}
\examples{lagcor(rnorm(10), rnorm(10), -3:3)}
