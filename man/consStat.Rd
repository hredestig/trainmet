\name{consStat}
\alias{consStat}
\title{Consensus statistic}
\usage{consStat(stat, members, cv=FALSE, fold=7)}
\description{Calculate a consensus statistic from a number of statistics that
each measure the association between a set of metabolites and a
set of transcripts.}
\details{With an matrix of association statistics between a set of
metabolites and transcripts this function uses a set of known
assocaition to perform a supervised summarization using the
PLS+CCA replacement for OPLS suggested by MacGregor. Briefly, a
PLS-DA model is calculated between using the known true
associations as guiding examples and PLS scores are
extracted. These are then rotated using  canonical covariate
analysis to obtain the component that is truly correlated to class
separation. The resulting score vector is the consensus statistic
(which is high for strong associations but has not other
interpretation).}
\value{A list with the consesus statistic (\code{cct}), the
correlation loadings (\code{ccp}) and possibly the cv predicted
statistic (\code{ccthat}).}
\references{Yu, H and MacGregor, JF. (2004) Post processing
methods (PLS-CCA): simple alternatives to pre-processing methods
(OSC-PLS). Chemometrics and Intelligent Laboratory Systems 73,
199-205.}
\author{Henning Redestig}
\arguments{\item{stat}{a matrix of association statistics (rows are transcripts)}
\item{members}{a logical vector indiating if the transcript is
associated to the metabolite(s).}
\item{cv}{Do cross-validation and return the cross-predicted}
\item{fold}{Number of CV-folds}}
\examples{stat <- matrix(rnorm(100*5), ncol=5)
members <- c(rep(1, 21), rep(0, 79))
consStat(stat, members)}
