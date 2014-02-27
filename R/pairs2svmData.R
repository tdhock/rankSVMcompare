### We do several preprocessing steps on the input data: \enumerate{
### \item centering and scaling. \item for all pairs i such that
### yi=-1, we flip Xi and Xip and set yi=1. \item for all pairs i such
### that yi=0, we generate new features Xi <- Xip and Xip <- Xi with
### corresponding new yi=0 labels. \item We take the difference Xip-Xi
### of the resulting scaled, flipped, augmented feature
### matrices. \item We map 0 -> -1 in the resulting label vector,
### creating an integer vector with elements in c(-1,1).}
pairs2svmData <- structure(function(Pairs){
  check.pairs(Pairs)
  ## First make scaled input features.
  scaled <- with(Pairs, scale(rbind(Xi, Xip)))
  sigma <- attr(scaled, "scaled:scale")
  include <- sigma != 0
  sigma <- sigma[include]
  if(length(sigma) == 0){
    stop("variance=0 for all input features")
  }
  mu <- attr(scaled, "scaled:center")[include]
  Xi <- scale(Pairs$Xi[,include,drop=FALSE], mu, sigma)
  Xip <- scale(Pairs$Xip[,include,drop=FALSE], mu, sigma)
  ## Then we put all nonzero y_i pairs on the y_i=1 side.
  yi <- Pairs$yi
  flip <- yi == -1
  tmp <- Xip[flip,]
  Xip[flip,] <- Xi[flip,]
  Xi[flip,] <- tmp
  yi[flip] <- 1
  ## And we duplicate the y_i=0 pairs, on the negative side.
  zero <- which(yi == 0)
  Xi <- rbind(Xi, Xip[zero,,drop=FALSE])
  Xip <- rbind(Xip, Xi[zero,,drop=FALSE])
  yi <- c(yi, rep(0, length(zero)))
  yi <- ifelse(yi==0, -1L, 1L)
  ##value<< data suitable for plugging into an SVM solver:
  list(center=mu, ##<< center of the input features.
       scale=sigma, ##<< scale of input features.
       Xi=Xi, ##<< inputs feature matrix.
       Xip=Xip, ##<< inputs feature matrix.
       yi=yi) ##<< outputs -1 -> 1, 0 -> -1, 1 -> 1.
  ##end<<
},ex=function(){
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
})
