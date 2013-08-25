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
  ## First make scaled input features Zi, Zip.
  scaled <- with(Pairs, scale(rbind(Xi, Xip)))
  mu <- attr(scaled, "scaled:center")
  sigma <- attr(scaled, "scaled:scale")
  Zi <- scale(Pairs$Xi, mu, sigma)
  Zip <- scale(Pairs$Xip, mu, sigma)
  Di <- Zip - Zi
  ## Then we put all nonzero y_i pairs on the y_i=1 side.
  yi <- Pairs$yi
  is.zero <- yi == 0
  Di[!is.zero,] <- Di[!is.zero,] * yi[!is.zero]
  yi[!is.zero] <- 1
  ## And we duplicate the y_i=0 pairs, on the negative side.
  Di.other <- Di
  Di.other[is.zero,] <- -Di.other[is.zero,]
  X <- rbind(Di,Di.other[is.zero,])
  yi.both <- c(yi, rep(0, sum(is.zero)))
  svm.y <- ifelse(yi.both==0, -1, 1)
  ##value<< data suitable for plugging into an SVM solver:
  list(scale=sigma, ##<< scale of input features.
       features=X, ##<< inputs: feature difference matrix Xip-Xi.
       labels=svm.y) ##<< outputs: -1 -> 1, 0 -> -1, 1 -> 1.
  ##end<<
},ex=function(){
  p <- list(Xi=rbind(3,0,1),
            Xip=rbind(0,-2,0),
            yi=as.integer(c(-1,1,0)))
  result <- pairs2svmData(p)
  ## Inequality pairs such that yi=1 or -1 are mapped to 1, and
  ## equality pairs such that yi=0 are duplicated and mapped to -1.
  stopifnot(result$labels == c(1,1,-1,-1))
  ## The duplicate equality features are multiplied by -1.
  stopifnot(result$features[3]*-1 == result$features[4])
})
