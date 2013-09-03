softCompareQP <- structure(function
### Fit a sparse kernel soft margin comparison model by using
### \code{\link{ksvm}} to solve the SVM dual problem. We first
### normalize the data using \code{\link{pairs2svmData}}, resulting in
### scaled n x p feature matrices Xi and Xip, with a new vector of n
### comparisons yi in c(-1,1). We then make the 2n x p matrix
### X=rbind(Xi,Xip) and calculate its 2n x 2n kernel matrix K. We then
### use \code{ksvm(M %*% K %*% t(M), yi)}, where the n x 2n matrix
### M=[-In In], and In is the n x n identity matrix. The 2n primal
### kernel coefficients are \code{a = t(M) %*% (yi*v)} where the n
### dual variables v are obtained from the ksvm solver. For a scaled
### p-vector x, the learned binary classification function is
### \eqn{f(x)=b+\sum_{i=1}^n a_i K(x_i, x) + a_{n+i} K(x_i', x)} and
### the learned ranking function is \eqn{r(x)=\sum_{i=1}^n -a_i/b
### K(x_i, x) - a_{n+i}/b K(x_i', x)}, where the scalar intercept/bias
### b is obtained from the ksvm solver.
(Pairs,
### see \code{\link{check.pairs}}.
 kernel=rbfdot(sigma=1),
### Kernel function, see \code{\link{dots}}.
 ...
### Passed to \code{\link{ksvm}}.
 ){
  if(is.character(kernel)){
    kernel <- get(kernel)
  }
  stopifnot(is.function(kernel))
  is.kernel <- function(k)inherits(k,"kernel")
  if(!is.kernel(kernel)){
    kernel <- kernel()
  }
  stopifnot(is.kernel(kernel))
  res <- pairs2svmData(Pairs)
  res$kernel <- kernel
  X.all <- with(res, rbind(Xi, Xip))
  K2n <- kernelMatrix(kernel, X.all)
  P <- ncol(X.all)
  N <- nrow(res$Xi)
  M <- rbind(diag(rep(-1,N)),diag(rep(1,N)))
  Kn <- as.kernelMatrix(t(M) %*% K2n %*% M)
  fit <- ksvm(Kn, res$yi, type="C-svc", ...)
  res$ksvm <- fit
  dual.var.times.y <- rep(0, N)
  dual.var.times.y[fit@SVindex] <- fit@coef[[1]]
  res$primal <- M %*% dual.var.times.y
  is.sv <- res$primal!=0
  res$sv <- list(X=X.all[is.sv,],
                 a=res$primal[is.sv])
  res$check <- function(X){
    stopifnot(ncol(X)==P)
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))
  }
  res$svm.f <- function(X){
    res$check(X)
    kernelMult(kernel, X, res$sv$X, res$sv$a)-fit@b
  }
  res$rank.scaled <- function(X){
    res$check(X)
    kernelMult(kernel, X, res$sv$X, res$sv$a/fit@b)
  }
  res$rank <- function(X){
    res$check(X)
    X.sc <- scale(X, res$center, res$scale)
    res$rank.scaled(X.sc)
  }
  res$predict <- function(Xi, Xip){
    rank.diff <- res$rank(Xip)-res$rank(Xi)
    ifelse(rank.diff < -1, -1L,
           ifelse(rank.diff > 1, 1L, 0L))
  }
  res
### Comparison model fit. You can do fit$rank(X) to get m numeric
### ranks for the rows of the unscaled m x p numeric matrix X. For two
### feature vectors xi and xip, we predict no significant difference
### if their absolute rank difference is less than 1. You can do
### fit$predict(Xi,Xip) to get m predicted comparisons in c(-1,0,1),
### for m by p numeric matrices Xi and Xip. Also, fit$scale and
### fit$center are the centers and scales of the input features, and
### fit$sv are the support vectors (in the scaled space). For a scaled
### matrix X, fit$svm.f(X) gives the values of the learned binary
### classification function f(x), and fit$rank.scaled(X) gives the
### values of the ranking function r(x).
},ex=function(){
  data(separable,package="rankSVMcompare")
  ## Add some noise to create a data set which is not separable.
  not <- separable
  set.seed(1)
  for(name in c("Xi","Xip")){
    not[[name]][,"distance"] <-
      not[[name]][,"distance"]+rnorm(nrow(not$Xi),sd=50)
  }
  arrow.df <- with(not, data.frame(Xi,Xip,yi))
  library(ggplot2)
  library(grid)
  ## This is the pattern in the training data.
  arrowPlot <- ggplot(,aes(distance, angle))+
    geom_segment(aes(xend=distance.1, yend=angle.1),
                 data=subset(arrow.df,yi==0))+
    geom_segment(aes(xend=distance.1, yend=angle.1),
                 data=subset(arrow.df,yi!=0),arrow=arrow())+
    facet_grid(.~yi)+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))
  print(arrowPlot)
  g.size <- 20
  d <- with(arrow.df, range(c(distance, distance.1)))
  a <- with(arrow.df, range(c(angle, angle.1)))
  X.grid <- as.matrix(
    expand.grid(distance=seq(d[1], d[2], l=g.size),
                angle=seq(a[1], a[2], l=g.size))
              )
  ## Fit some soft-margin linear comparison models.
  unscaledGrid <- data.frame()
  scaledGrid <- data.frame()
  seg.df <- data.frame()
  sv.df <- data.frame()
  for(cost in c(0.01, 1)){
    fit <- softCompareQP(not, kernel="vanilladot", C=cost)
    point.df <- with(fit, data.frame(Xip-Xi, yi))
    d <- range(point.df$distance)
    a <- range(point.df$angle)
    sc.grid <- as.matrix(expand.grid(distance=seq(d[1],d[2],l=g.size),
                                     angle=seq(a[1],a[2],l=g.size)))
    for(fun.name in c("svm.f","rank.scaled")){
      fun <- fit[[fun.name]]
      value <- fun(sc.grid)
      scaledGrid <- rbind(scaledGrid,{
        data.frame(cost, sc.grid, value, fun.name)
      })
    }
    dual.var <- fit$ksvm@coef[[1]]
    support.vectors <- with(fit, {
      data.frame((Xip-Xi)[ksvm@SVindex,], cost,
                 sv.type=ifelse(abs(dual.var)==max(dual.var),"slack","margin"))
    })
    sv.df <- rbind(sv.df, support.vectors)
    unscaledGrid <- rbind(unscaledGrid, {
      data.frame(X.grid, rank=fit$rank(X.grid), cost)
    })
  }
  library(directlabels)
  arrowContour <- arrowPlot+
    geom_contour(aes(z=rank, colour=..level..), size=1.5, data=unscaledGrid)+
    geom_dl(aes(z=rank, colour=..level.., label=..level..), 
            data=unscaledGrid, method="bottom.pieces", stat="contour")+
    facet_grid(cost~yi)
  print(arrowContour)
  ## Since we learned a linear comparison model we can also interpret
  ## the model in the difference space.
  brks <- c(-2,-1,0,1,2)
  diffPlot <- ggplot()+
    geom_point(aes(distance, angle, colour=factor(yi)), data=point.df)+
    geom_contour(aes(distance, angle, z=value), size=1.5,
                 data=scaledGrid, colour="grey", breaks=brks)+
    geom_dl(aes(distance, angle, z=value, label=..level..), colour="grey",
            data=scaledGrid, method="bottom.pieces",
            breaks=brks, stat="contour")+
    ## geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
    ##                  linetype=line),data=seg.df)+
    geom_point(aes(distance, angle, shape=sv.type), data=sv.df,size=5)+
    scale_linetype_manual(values=c(margin="dashed",decision="solid"))+
    scale_shape_manual(values=c(margin=13,slack=3))+
    facet_grid(fun.name~cost,labeller=function(var, val){
      if(var=="cost"){
        sprintf("C = %s",as.character(val))
      }else as.character(val)
    })+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))+
    ggtitle("Support vector comparison model (grey level lines)")
  print(diffPlot)
  ## Fit soft-margin non-linear comparison models for some cost and
  ## kernel width parameters.
  fold <- sample(rep(1:4, l=nrow(not$Xi)))
  set.ids <- list(train=1:2, validation=3, test=4)
  sets <- list()
  for(set.name in names(set.ids)){
    i <- fold %in% set.ids[[set.name]]
    sets[[set.name]] <- with(not, list(Xi=Xi[i,], Xip=Xip[i,], yi=yi[i]))
  }
  grid.df <- data.frame()
  err.df <- data.frame()
  model.i <- 1
  fits <- list()
  for(cost in c(0.1, 10)){
    for(sigma in c(1/4, 4)){
      fit <- softCompareQP(sets$train, kernel=rbfdot(sigma), C=cost)
      fits[[model.i]] <- fit
      yhat.vali <- with(sets$validation, fit$predict(Xi, Xip))
      ## Error on the validation set used for model selection:
      table(yhat.vali, sets$validation$yi)
      err <- FpFnInv(yhat.vali, sets$validation$yi)
      err.df <- rbind(err.df, data.frame(err, cost, sigma, model.i))
      f <- fit$rank(X.grid)
      grid.df <- rbind(grid.df, data.frame(X.grid, f, cost, sigma, model.i))
      model.i <- model.i+1
    }
  }
  nonlinear <- arrowPlot+
    geom_contour(aes(distance, angle, z=f), size=1.5,
                 data=grid.df, colour="grey")+
    geom_dl(aes(distance, angle, z=f, label=..level..), colour="grey",
            data=grid.df, method="bottom.pieces", stat="contour")+
    facet_grid(model.i~yi,labeller=function(var, val){
      sprintf("%s = %s",var,as.character(val))
    })+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))+
    ggtitle("Support vector comparison model (grey level curves)")
  print(nonlinear)
  ## Model selection:
  fit <- fits[[which.min(err.df$error)]]
  ## Out of sample prediction error:
  test.yhat <- with(sets$test, fit$predict(Xi,Xip))
  with(sets$test, table(labels=yi, predictions=test.yhat))
  with(sets$test, FpFnInv(yi, test.yhat))
  ## With this setup, we have no prediction error!

  ## graph data.
  if(FALSE){
    data(graphs)
    set.ids <- list(train=1:3,
                    validation=4:5,
                    test=6:7)
    sets <- list()
    for(set.name in names(set.ids)){
      i <- set.ids[[set.name]]
      L <- graphs[i]
      sets[[set.name]] <-
        list(Xi=do.call(rbind, lapply(L, "[[", "Xi")),
             Xip=do.call(rbind, lapply(L, "[[", "Xip")),
             yi=do.call(c, lapply(L, "[[", "yi")))
    }
    err.df <- data.frame()
    model.i <- 1
    fits <- list()
    for(cost in 2^( (-5):5 )){
      for(sigma in 2^( (-5):5 )){
        cat(sprintf("model %4d sigma=%10f cost=%10f\n", model.i, sigma, cost))
        fit <- softCompareQP(sets$train, kernel=rbfdot(sigma=sigma), C=cost)
        fits[[model.i]] <- fit
        yhat.vali <- with(sets$validation, fit$predict(Xi, Xip))
        ## Error on the validation set used for model selection:
        table(yhat.vali, sets$validation$yi)
        err <- FpFnInv(yhat.vali, sets$validation$yi)
        err.df <- rbind(err.df, data.frame(err, cost, sigma, model.i))
        model.i <- model.i+1
      }
    }
    ## Table and heatmap of validation error.
    err.mat <- sapply(split(err.df, err.df$cost), "[[", "error")
    print(err.mat)
    heat <- ggplot(err.df, aes(log2(cost), log2(sigma)))+
      geom_tile(aes(fill=error))
    print(heat)
    ## Model selection and test error.
    best.i <- err.df[which.min(err.df$error),"model.i"]
    fit <- fits[[best.i]]
    test.pred <- with(sets$test, fit$predict(Xi, Xip))
    table(labels=sets$test$yi, predictions=test.pred)
    stats <- FpFnInv(sets$test$yi, test.pred)
    with(stats, error/count) #proportion wrong = test error
  }
})
