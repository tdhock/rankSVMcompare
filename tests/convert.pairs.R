library(rankSVMcompare)

p <- list(Xi=cbind(var=c(3,0,1)),
          Xip=cbind(var=c(0,-2,0)),
          yi=as.integer(c(-1,1,0)))
result <- pairs2svmData(p)
## Inequality pairs such that yi=1 or -1 are mapped to 1, and
## equality pairs such that yi=0 are duplicated and mapped to -1.
stopifnot(result$labels == c(1,1,-1,-1))
## The duplicate equality features are flipped.
stopifnot(result$Xi[3]==result$Xip[4])
stopifnot(result$Xi[4]==result$Xip[3])

for(i in seq_along(p)){
  colnames(p[[i]]) <- NULL
}
result <- tryCatch({
  pairs2svmData(p)
}, error=function(e){
  e$mes
})
stopifnot(is.character(result))
stopifnot(length(result) == 1)
stopifnot(result == "input feature matrix columns must have names")
