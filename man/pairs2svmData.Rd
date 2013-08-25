\name{pairs2svmData}
\alias{pairs2svmData}
\title{pairs2svmData}
\description{We do several preprocessing steps on the input data: \enumerate{
\item centering and scaling. \item for all pairs i such that
yi=-1, we flip Xi and Xip and set yi=1. \item for all pairs i such
that yi=0, we generate new features Xi <- Xip and Xip <- Xi with
corresponding new yi=0 labels. \item We take the difference Xip-Xi
of the resulting scaled, flipped, augmented feature
matrices. \item We map 0 -> -1 in the resulting label vector,
creating an integer vector with elements in c(-1,1).}}
\usage{pairs2svmData(Pairs)}
\arguments{
  \item{Pairs}{
}
}

\value{data suitable for plugging into an SVM solver:
\item{scale}{scale of input features.}
\item{features}{inputs: feature difference matrix Xip-Xi.}
\item{labels}{outputs: -1 -> 1, 0 -> -1, 1 -> 1.}}

\author{Toby Dylan Hocking}





