\name{makeMetLabels2}
\alias{makeMetLabels2}
\title{Make unique non-mising labels for a set of metabolites}
\usage{
  makeMetLabels2(object,
    targets = c("preferred", "synonym", "lineno"))
}
\arguments{
  \item{object}{an ExpressionSet}

  \item{targets}{the selected columns in fData to search
  in.}
}
\value{
  character vector with non-missing unique elements
}
\description{
  Labels are made by taking the target columns from fData
  and filling labels sequentially from the columns in fData
  prioritizing the selected columns. Labels are then made
  unique and any remaining missing values are replaced by a
  constructed identifier.
}
\examples{
data(mdat)
makeMetLabels2(mdat)
}
\author{
  henning
}

