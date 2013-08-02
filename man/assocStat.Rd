\name{assocStat}
\alias{assocStat}
\title{Other statistics}
\usage{assocStat(tdat, mdat, method=c("pearson", "spearman", "lagpearson",
    "hmm"), absolute=TRUE, ...)}
\description{Compute other (than HMM based) statistics for associations between
transcripts and metabolites.}
\value{A matrix with statistics}
\author{Henning Redestig}
\arguments{\item{tdat}{transcripts ExpressionSet}
\item{mdat}{metabolites ExpressionSet}
\item{method}{the method to use}
\item{absolute}{logical indicting if inverted patterns also be
considered. For ordinary Pearson, this implies concatenating the
absolute correlations to the result.}
\item{...}{arguments passed on to the chosen method, i.e.,
or \code{lagcor} or \code{trainrank} (for the HMM)}}
\examples{data(tdat)
data(mdat)
library(org.At.tair.db)
stat <- assocStat(tdat, mdat, use="pairwise")
members <- rep(FALSE, nrow(tdat))
members[match(get("00052", org.At.tairPATH2TAIR), fData(tdat)$agi, nomatch=0)] <- TRUE
cstat <- consStat(stat, members)
boxplot(cstat$cct~members)}
