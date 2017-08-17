\name{FpFnInv}
\alias{FpFnInv}
\title{FpFnInv}
\description{Count false positive, false negative, and inversions. Each of
these is considered an error. Both the real labels and the
predictions should be in c(-1,0,1).}
\usage{FpFnInv(true, pred)}
\arguments{
  \item{true}{Actual labels.}
  \item{pred}{Model predictions.}
}

\value{data.frame with 1 row and 5 columns: false.positive,
false.negative, inversion, error, count.}

\author{Toby Dylan Hocking}




\examples{

values <- c(-1,0,1)
x <- expand.grid(true=values, pred=values)
(err <- FpFnInv(x$true, x$pred))
stopifnot(err$err == 6)
stopifnot(err$false.positive == 2)
stopifnot(err$false.negative == 2)
stopifnot(err$inversion == 2)
stopifnot(err$count == nrow(x))

}
