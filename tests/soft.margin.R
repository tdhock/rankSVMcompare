library(rankSVMcompare)

data(separable)

no.var <- cbind(zero=rep(0, nrow(separable$Xi)),
                one=rep(1, nrow(separable$Xi)))
bad.pairs <- list(Xi=no.var, Xip=no.var, yi=separable$yi)
good.bad <- with(separable, {
  list(Xi=cbind(Xi, no.var),
       Xip=cbind(Xip, no.var),
       yi=yi)
})
## Test data that does not have all the required columns.
not.enough <- with(separable, {
  list(Xi=cbind(Xi[,1,drop=FALSE], no.var),
       Xip=cbind(Xip[,1,drop=FALSE], no.var),
       yi=yi)
})
## Test data that does not have named columns.
not.named <- good.bad
for(i in seq_along(not.named)){
  colnames(not.named[[i]]) <- NULL
}
check <- function(Pairs, error=NULL){
  list(data=Pairs, error=error)
}
## Check that these data sets can be trained, or result in an error.
to.check <-
  list(check(good.bad),
       check(separable),
       check(not.named, error="input feature matrix columns must have names"),
       check(bad.pairs, error="variance=0 for all input features"))
## For each data set that can be trained, check that these data sets
## can be tested, or result in an error.
to.test <-
  list(check(good.bad),
       check(separable),
       check(bad.pairs, "missing features distance, angle"),
       check(not.enough, "missing features angle"),
       check(not.named, "missing features distance, angle"))
for(L in to.check){
  fit.or.err <- tryCatch({
    softCompareQP(L$data)
  }, error=function(e){
    e$mes
  })
  if(is.null(L$error)){ ## result.
    stopifnot(is.list(fit.or.err))
    pred <- fit.or.err$predict
    stopifnot(is.function(pred))
    for(test in to.test){
      labels.or.err <- tryCatch({
        with(test$data, pred(Xi, Xip))
      }, error=function(e){
        e$mes
      })
      if(is.null(test$error)){ #result.
        stopifnot(labels.or.err %in% c(-1, 0, 1))
        stopifnot(length(labels.or.err) == nrow(test$data))
      }else{ #error.
        if(!identical(test$error, labels.or.err)){
          print(test$error)
          print(labels.or.err)
          stop("expected first, but got second")
        }
      }
    }
  }else{ ##error.
    if(!identical(L$error, fit.or.err)){
      print(L$error)
      print(fit.or.err)
      stop("expected first, but got second")
    }
  }
}
