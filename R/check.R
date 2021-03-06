### Check for valid format of paired features and comparisons. The
### format is a list with n x p numeric matrices Xi and Xip, and yi is
### an integer vector of length n that takes values in c(-1,0,1):
### \itemize{ \item \eqn{y_i=-1} means that \eqn{x_i} is better than
### \eqn{x_i'}.\item \eqn{y_i=0} means that there is no significant
### difference between \eqn{x_i} and \eqn{x_i'}. \item \eqn{y_i=1}
### means that \eqn{x_i'} is better than \eqn{x_i}.}
check.pairs <- function(Pairs){
  N <- length(Pairs$yi)
  stopifnot(Pairs$yi %in% c(-1,0,1))
  P <- ncol(Pairs$Xi)
  for(v.name in c("Xi", "Xip")){
    mat <- Pairs[[v.name]]
    stopifnot(is.matrix(mat))
    stopifnot(is.numeric(mat))
    stopifnot(nrow(mat) == N)
    stopifnot(ncol(mat) == P)
    if(is.null(colnames(mat))){
      stop("input feature matrix columns must have names")
    }
  }
### Nothing.
}
