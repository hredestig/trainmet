\name{colorLineplot}
\alias{colorLineplot}
\title{Color line plot}
\usage{colorLineplot(tdat, mdat, stat, members, toprat=0.05, rows, cols,
    cex.text=1.1, ...)}
\description{Make a line plot of transcription and metabolite abundnance
trajectories over a time-course and color the genes after their
ranked co-response to the metabolite.}
\details{Can be used to visualize how highly each metabolite ranks the
different transcripts. Ranks are relative to the all genes but
only the selected genes are plotted.}
\value{nothing}
\author{Henning Redestig}
\arguments{\item{tdat}{the gene expression data. Typically the whole dataset.}
\item{mdat}{the metabolite data. Only the metabolites of
interest.}
\item{stat}{pre-calculated statistics that ranks all transcripts
after their association with the metabolites (genes in rows)}
\item{members}{numerical indicating the genes that should be plotted}
\item{toprat}{the top ratio of genes that should receive the
highest color intensity}
\item{rows}{display layout, rows}
\item{cols}{display layout,k columns}
\item{cex.text}{expansion factor for the text}
\item{...}{passed on to \code{par}}}
\examples{data(tdat)
data(mdat)
library(org.At.tair.db)
pearson <- assocStat(tdat, mdat, absolute=FALSE)
members <- na.omit(match(get("00052", org.At.tairPATH2TAIR), fData(tdat)$agi))
colorLineplot(tdat, mdat, pearson, members, rows=2, cols=3)}
