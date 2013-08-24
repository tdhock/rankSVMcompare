### Check for valid format of paired features and comparisons.
check.pairs <- function(Pairs){
  stopifnot(is.integer(Pairs$yi))
  stopifnot(is.numeric(Pairs$Xi))
  stopifnot(is.numeric(Pairs$Xip))
  N <- length(Pairs$yi)
  stopifnot(Pairs$yi %in% c(-1,0,1))
  stopifnot(nrow(Pairs$Xi) == N)
  stopifnot(nrow(Pairs$Xip) == N)
  P <- ncol(Pairs$Xi)
  stopifnot(P == ncol(Pairs$Xip))
}
