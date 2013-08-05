\name{trainrank}
\alias{trainrank}
\title{The metabolites/transcripts HMM association program}
\usage{
  trainrank(mdat, tdat, states = 3, cyclic = FALSE,
    LD_LIBRARY_PATH = "/usr/local/lib/",
    trainrank = system.file("extdata", "train_rank.py", package = "trainmet"),
    python = "python -W ignore", verbose = FALSE,
    cleanOnExit = TRUE, scale = TRUE)
}
\arguments{
  \item{mdat}{matrix with the metabolite data (time points
  in columns)}

  \item{tdat}{matrix with the transcript data (time points
  in columns)}

  \item{states}{number of states}

  \item{cyclic}{use cyclic model}

  \item{LD_LIBRARY_PATH}{path to ghmm shared object files}

  \item{trainrank}{where to find the train_rank.py script}

  \item{python}{the command use to call python}

  \item{verbose}{be verbose}

  \item{cleanOnExit}{if true, remove temp files and created
  directories on exit}

  \item{scale}{logical indicating if the data should be
  scaled to unitvariance and centered. HMMs are very
  sensitive to this and scaling is strongly recommended.}
}
\value{
  a list with results for negative and positive modes
}
\description{
  R wrapper for the train_rank.py script.
}
\details{
  For larger numbers of metabolites (>50) calling this
  function will take a while.
}
\examples{
trainrank(rbind(sample(1:10)), rbind(rnorm(10)))
metabolite <- jitter(c(1,1,1,2,2,2,6,6,6,6,6,2,2,2))
associated_gene <- jitter(c(2,2,2,6,6,6,6,6,2,2,2,3,3,3))
unassociated_genes <- t(replicate(50, rnorm(14)))
mdat <- rbind(metabolite)
tdat <- rbind(associated_gene, unassociated_genes)
rownames(mdat) <- rownames(tdat) <- NULL
matplot(cbind(metabolite, associated_gene), type='b')
tt <- trainrank(mdat, tdat)
barplot(tt$positive[,1])
tt$plag[1,1]
}
\author{
  Henning Redestig
}
\references{
  Redestig, H. and Costa, IG (submitted)
}

