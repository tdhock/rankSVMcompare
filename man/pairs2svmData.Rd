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
\item{center}{center of the input features.}
\item{scale}{scale of input features.}
\item{Xi}{inputs feature matrix.}
\item{Xip}{inputs feature matrix.}
\item{yi}{outputs -1 -> 1, 0 -> -1, 1 -> 1.}}

\author{Toby Dylan Hocking}




\examples{

p <- list(Xi=cbind(var=c(3,0,1)),
          Xip=cbind(var=c(0,-2,0)),
          yi=as.integer(c(-1,1,0)))
(result <- pairs2svmData(p))
## Inequality pairs such that yi=1 or -1 are mapped to 1, and
## equality pairs such that yi=0 are duplicated and mapped to -1.
stopifnot(result$labels == c(1,1,-1,-1))
## The duplicate equality features are flipped.
stopifnot(result$Xi[3]==result$Xip[4])
stopifnot(result$Xi[4]==result$Xip[3])

}
