\name{guessComp2}
\alias{guessComp2}
\title{Guess optimal number of components...}
\usage{guessComp2(fit, y, lev=0.05, minq2=0)}
\description{Guess optimal number of components}
\details{Use some heuristics to guess the best number of components for PLS
model. Specifically, take the best component yielding the best Q2
and then remove components as long as there is no sudden decrease
in Q2.}
\value{a list with R2, Q2 and the recommended number of components}
\author{Henning Redestig}
\arguments{\item{fit}{a \code{plsr} model}
\item{y}{the y data}
\item{lev}{the ratio of Q2 that is considered small a decrease.}
\item{minq2}{the minimum Q2 to consider, less than this and zero is
returned instead.}}
\examples{data(oliveoil)
sens.pls <- plsr(sensory~chemical, validation="CV", ncomp=5, data=oliveoil)
with(oliveoil, guessComp2(sens.pls, sensory))}
