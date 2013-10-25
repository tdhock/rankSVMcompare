FpFnInv <- structure(function
### Count false positive, false negative, and inversions. Each of
### these is considered an error. Both the real labels and the
### predictions should be in c(-1,0,1).
(true,
### Actual labels.
 pred
### Model predictions.
 ){
  stopifnot(length(true)==length(pred))
  stopifnot(c(true, pred) %in% c(-1,0,1))
  fp <- sum(true == 0 & pred != 0)
  inv <- sum(true %in% c(-1,1) & pred != true & pred %in% c(-1,1))
  fn <- sum(true !=0 & pred==0)
  err <- sum(true != pred)
  stopifnot(err == fp+inv+fn)
  data.frame(false.positive=fp,
             false.negative=fn,
             inversion=inv,
             error=err,
             equality=sum(true==0),
             inequality=sum(true!=0),
             count=length(true))
### data.frame with 1 row and 5 columns: false.positive,
### false.negative, inversion, error, count.
},ex=function(){
  values <- c(-1,0,1)
  x <- expand.grid(true=values, pred=values)
  err <- FpFnInv(x$true, x$pred)
  stopifnot(err$err == 6)
  stopifnot(err$false.positive == 2)
  stopifnot(err$false.negative == 2)
  stopifnot(err$inversion == 2)
  stopifnot(err$count == nrow(x))
})
