softCompareQP <- structure(function
### Fit a soft-margin comparison model by using \code{\link{ksvm}}
### (libsvm) to solve the dual SVM problem. TODO: explain optimization
### problem.
(Pairs,
### see \code{\link{check.pairs}}.
 ...
### Passed to \code{\link{ksvm}}.
 ){
  res <- pairs2svmData(Pairs)
  X <- res$features
  P <- ncol(X)
  y <- res$labels
  diff.df <- data.frame(X, y)
  fit <- ksvm(X, y, type="C-svc", scaled=FALSE, ...)
  res$ksvm <- fit
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
    f2 <- predict(fit, -X.sc, type="decision")
    ##TODO: FIX! This works for the linear but nonlinear case!
    (f.svm-f2)/res$gamma
  }
  ## Find the optimal nonlinear scaling, via grid search?
  ## thresh <- function(x)ifelse(x < -1, 1L, ifelse(x > 1, 1L, -1L))
  ##table(thresh(res$rank(X)), sign(predict(fit, X, type="decision")))
  res$gamma <- fit@b*2
  res$predict <- function(Xi, Xip){
    for(X in list(Xi, Xip)){
      res$check(X)
    }
    D <- Xip-Xi
    ## TODO: fix, check what makes sense here for the nonlinear case!
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
  g.size <- 20
  X.grid <- as.matrix(with(point.df, {
    expand.grid(distance=seq(min(distance), max(distance), l=g.size),
                angle=seq(min(angle), max(angle), l=g.size))
  }))
  ## Fit 3 soft-margin linear comparison models.
  grid.df <- data.frame()
  seg.df <- data.frame()
  sv.df <- data.frame()
  for(cost in c(0.01, 1, 100)){
    fit <- softCompareQP(not, kernel="vanilladot", C=cost)
    mu <- fit$margin
    w <- fit$weight
    arange <- range(point.df$angle)
    seg <- function(v, line){
      d <- (v-w[2]*arange)/w[1]
      data.frame(t(c(distance=d, angle=arange)), line, cost)
    }
    seg.df <- rbind(seg.df,
                    seg(1-mu,"margin"),
                    seg(1+mu,"margin"),
                    seg(-1-mu,"margin"),
                    seg(-1+mu,"margin"),
                    seg(1,"decision"),
                    seg(-1,"decision"))
    support.vectors <- with(fit, {
      data.frame(t(t(sv$X)*fit$scale), cost,
                 sv.type=ifelse(abs(sv$a)==max(sv$a),"slack","margin"))
    })
    sv.df <- rbind(sv.df, support.vectors)
    f <- fit$rank(X.grid)
    grid.df <- rbind(grid.df, data.frame(X.grid, f, cost))
  }
  library(directlabels)
  library(grid)
  pmodel <- p+
    geom_contour(aes(distance, angle, z=f), size=1.5,
                 data=grid.df, colour="grey")+
    geom_dl(aes(distance, angle, z=f, label=..level..), colour="grey",
            data=grid.df, method="bottom.pieces", stat="contour")+
    geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                     linetype=line),data=seg.df)+
    geom_point(aes(distance, angle, shape=sv.type), data=sv.df)+
    scale_linetype_manual(values=c(margin="dashed",decision="solid"))+
    scale_shape_manual(values=c(margin=13,slack=3))+
    facet_grid(.~cost,labeller=function(var, val){
      sprintf("C = %s",as.character(val))
    })+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))+
    ggtitle("Support vector comparison model (black and grey level lines)")
  print(pmodel)
  ## Fit soft-margin non-linear comparison models for a 2 x 2 grid of
  ## cost and kernel width parameters.
  grid.df <- data.frame()
  sv.df <- data.frame()
  for(cost in c(0.01, 1, 100)){
    for(sigma in c(1/2, 2)){
      fit <- softCompareQP(not, kernel="rbfdot", C=cost, kpar=list(sigma=sigma))
      support.vectors <- with(fit, {
        data.frame(t(t(sv$X)*fit$scale), cost, sigma,
                   sv.type=ifelse(abs(sv$a)==max(sv$a),"slack","margin"))
      })
      sv.df <- rbind(sv.df, support.vectors)
      f <- fit$rank(X.grid)
      grid.df <- rbind(grid.df, data.frame(X.grid, f, cost, sigma))
    }
  }
  library(directlabels)
  library(grid)
  nonlinear <- p+
    geom_contour(aes(distance, angle, z=f), size=1.5,
                 data=grid.df, colour="grey")+
    geom_dl(aes(distance, angle, z=f, label=..level..), colour="grey",
            data=grid.df, method="bottom.pieces", stat="contour")+
    geom_point(aes(distance, angle, shape=sv.type), data=sv.df)+
    scale_linetype_manual(values=c(margin="dashed",decision="solid"))+
    scale_shape_manual(values=c(margin=13,slack=3))+
    facet_grid(sigma~cost,labeller=function(var, val){
      sprintf("%s = %s",var,as.character(val))
    })+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))+
    ggtitle("Support vector comparison model (grey level curves)")
  print(nonlinear)
})
