\name{separable}
\alias{separable}
\docType{data}
\title{
Separable data for hard margin comparison models
}
\description{
  Two features (angle and distance) were measured for several layouts of
  a graph with 2 nodes. Each pair of layouts is represented by a row of
  Xi and Xip, and the corresponding element of the labels yi.
  The goal is analyzing
  these data is to find a function that accurately predicts the label yi
  given a pair of features Xi and Xip.
}
\usage{data(separable)}
\format{
  List of n=189 by p=2 numeric feature matrices Xi and Xip, with
  n corresponding comparisons yi. See \code{\link{check.pairs}}.
}
\keyword{datasets}
