\name{getcvseg}
\alias{getcvseg}
\title{Class balanced CV segments}
\usage{
  getcvseg(x, fold = 7, seed = 123)
}
\arguments{
  \item{x}{a factor, character or numeric vector that
  describes class membership of a set of items.}

  \item{fold}{the desired number of segments}

  \item{seed}{randomization seed for reproducibility}
}
\value{
  a matrix where each row is CV segment
}
\description{
  Get cross-validation segments that have (as far as
  possible) the same ratio of all classes.
}
\examples{
seg <- getcvseg(iris$Species, 10)
apply(seg, 1, function(s) table(iris$Species[s]))
}
\author{
  Henning Redestig
}

