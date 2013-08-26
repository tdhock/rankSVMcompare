softCompareQP <- structure(function
### Fit a soft-margin comparison model by using kernlab::ksvm (libsvm)
### to solve the dual SVM problem. TODO: explain optimization problem.
(Pairs,
### see \code{\link{check.pairs}}.
 ...
### Passed to ksvm.
 ){
  res <- pairs2svmData(Pairs)
  X <- res$features
  P <- ncol(X)
  y <- res$labels
  diff.df <- data.frame(X, y)
  fit <- ksvm(X, y, type="C-svc", scaled=FALSE, ...)
  res$ksvm <- fit
  ## Make a grid and evaluate the learned fun on it.
  ## g.args <- list()
  ## for(name in colnames(X)){
  ##   x <- X[,name]
  ##   g.args[[name]] <- seq(min(x),max(x),l=20)
  ## }
  ## g <- do.call(expand.grid, g.args)
  ## g.unscaled <- t(t(g) * res$scale)
  ## f.svm <- predict(fit, g, type="decision")
  ## f.rank <- f.svm/fit@b+1
  ## res$grid.scaled <- data.frame(g, f.svm, f.rank)
  ## res$grid.unscaled <- data.frame(g.unscaled, f.svm, f.rank)
  res$sv <- list(X=fit@xmatrix[[1]],
                 a=fit@coef[[1]])
  res$weight.svm <- with(res$sv, colSums(X * a))
  res$margin <- 1/fit@b # not negative here, since f(x)=sum_i a_i x_i - b.
  res$weight <- with(res, margin*weight.svm/res$scale)
  res$check <- function(X){
    stopifnot(ncol(X)==P)
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))
  }
  res$rank <- function(X){
    res$check(X)
    X.sc <- t(t(X)/res$scale)
    f.svm <- predict(fit, X.sc, type="decision")
    1 + f.svm/fit@b
  }
  res$predict <- function(Xi, Xip){
    for(X in list(Xi, Xip)){
      res$check(X)
    }
    D <- Xip-Xi
    rank.diff <- res$rank(D)
    ifelse(rank.diff < -1, -1L,
           ifelse(rank.diff > 1, 1L, 0L))
  }
  res
},ex=function(){
  data(separable,package="rankSVMcompare")
  ## Add some noise to create a data set which is not separable.
  not <- separable
  set.seed(1)
  for(name in c("Xi","Xip")){
    not[[name]][,"distance"] <-
      not[[name]][,"distance"]+rnorm(nrow(not$Xi),sd=50)
  }
  point.df <- with(not, data.frame(Xip-Xi, yi))
  library(ggplot2)
  p <- ggplot()+
    geom_point(aes(distance, angle, colour=factor(yi)), data=point.df)
  print(p)
  ## Fit the soft-margin comparison model.
  fit <- softCompareQP(not, kernel="vanilladot")
  mu <- fit$margin
  w <- fit$weight
  arange <- range(point.df$angle)
  seg <- function(v, line){
    d <- (v-w[2]*arange)/w[1]
    data.frame(t(c(distance=d, angle=arange)), line)
  }
  seg.df <- rbind(seg(1-mu,"margin"),
                  seg(1+mu,"margin"),
                  seg(-1-mu,"margin"),
                  seg(-1+mu,"margin"),
                  seg(1,"decision"),
                  seg(-1,"decision"))
  sv.df <- with(fit, {
    data.frame(t(t(sv$X)*fit$scale),
               sv.type=ifelse(abs(sv$a)==max(sv$a),"error","margin"))
  })
  g.size <- 20
  X.grid <- with(point.df, {
    expand.grid(distance=seq(min(distance), max(distance), l=g.size),
                angle=seq(min(angle), max(angle), l=g.size))
  })
  X.grid$f <- fit$rank(as.matrix(X.grid))
  library(directlabels)
  p <- ggplot()+
    geom_contour(aes(distance, angle, z=f), size=1.5,
                 data=X.grid, colour="grey")+
    geom_dl(aes(distance, angle, z=f, label=..level..), colour="grey",
            data=X.grid, method="bottom.pieces", stat="contour")+
    ## geom_contour(aes(distance, angle, z=f.rank),
    ##              data=fit$grid.unscaled, colour="grey")+
    ## geom_dl(aes(distance, angle, z=f.rank, label=..level..), colour="grey",
    ##         data=fit$grid.unscaled, method="bottom.pieces", stat="contour")+
    geom_point(aes(distance, angle, colour=factor(yi)), data=point.df)+
    geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                     linetype=line),data=seg.df)+
    geom_point(aes(distance, angle, shape=sv.type), data=sv.df)+
    scale_linetype_manual(values=c(margin="dashed",decision="solid"))+
    scale_shape_manual(values=c(margin=13,error=3))
  print(p)
})
