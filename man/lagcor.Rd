\name{lagcor}
\alias{lagcor}
\title{Lagged pearson}
\usage{
  lagcor(x, y, lags = -2:2, ...)
}
\arguments{
  \item{x}{a numeric vector}

  \item{y}{a numeric vector}

  \item{lags}{lags to consider for the x-vector}

  \item{...}{passed on to cor}
}
\value{
  the best correlation over the considered lags
}
\description{
  (best) Lagged pearson
}
\examples{
lagcor(rnorm(10), rnorm(10), -3:3)
}
\author{
  henning
}

