hardCompareLP <- structure(function
### Fit a linear hard margin comparison model to linearly separable
### data. The linear program (LP) is max_{\eqn{\mu\in R}, \eqn{w\in
### R^p}} mu subject to these constraints: \itemize{ \item for all i
### such that \eqn{y_i=0}, \eqn{\mu \le 1-w'(x_i'-x_i)} and \eqn{\mu
### \le 1+w'(x_i'-x_i)}. \item for all i such that \eqn{y_i\ne 0},
### \eqn{\mu \le -1 + w'(x_i'-x_i)y_i}.}
(Pairs
### see \code{\link{check.pairs}}.
 ){
  check.pairs(Pairs)
  Di <- with(Pairs, Xip-Xi)
  N <- nrow(Di)
  P <- ncol(Di)
  vars <- make.ids(margin=1, weight=P)
  constraints <- list()
  for(i in 1:N){
    if(Pairs$yi[i] == 0){
      right.side <- -1
      yi.vec <- c(-1,1)
    }else{
      right.side <- 1
      yi.vec <- Pairs$yi[i]
    }
    di <- Di[i,]
    for(yi in yi.vec){
      const <- with(vars,{
        weight*di*yi + margin*-1 >= right.side
      })
      constraints <- c(constraints,list(const))
    }
  }
  n.vars <- length(unlist(vars))
  d <- rep(0, n.vars)
  d[vars$margin] <- 1
  result <- run.lpSolveAPI(vars, d, constraints)
  result$check <- function(X){
    stopifnot(ncol(X)==P)
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))
  }
  result$rank <- function(X){
    result$check(X)
    X %*% result$weight
  }
  result$predict <- function(Xi, Xip){
    for(X in list(Xi, Xip)){
      result$check(X)
    }
    D <- Xip-Xi
    rank.diff <- result$rank(D)
    ifelse(rank.diff < -1, -1L,
           ifelse(rank.diff > 1, 1L, 0L))
  }
  result
### Comparison model fit. You can do fit$rank(X) to get m numeric
### ranks for the rows of the m x p numeric matrix X. For two feature
### vectors xi and xip, we predict no significant difference if their
### absolute rank difference is less than 1. You can do
### fit$predict(Xi,Xip) to get m predicted comparisons in c(-1,0,1),
### for m by p numeric matrices Xi and Xip. Also, fit$weight is the
### optimal vector of p numeric weights, and if fit$margin is positive
### then the data are separable.
},ex=function(){
  data(separable)
  sol <- hardCompareLP(separable)
  ## check to make sure we have perfect prediction.
  y.hat <- with(separable, sol$predict(Xi, Xip))
  stopifnot(separable$yi == y.hat)
  ## This should also be the same:
  fxdiff <- with(separable, sol$rank(Xip)-sol$rank(Xi))
  y.hat2 <- ifelse(fxdiff < -1, -1L,
                   ifelse(fxdiff > 1, 1L, 0L))
  stopifnot(y.hat == y.hat2)

  ## Calculate which points are on the margin.
  margin <- ifelse(separable$yi==0,{
    1-abs(fxdiff)
  },{
    -1 + separable$yi * fxdiff
  })
  on.margin <- abs(margin - sol$margin)<1e-6
  diffs <- with(separable, {
    data.frame(Xip-Xi, yi,
               constraint=ifelse(on.margin, "active", "inactive"))
  })

  ## Calculate the decision and margin lines.
  arange <- range(diffs$angle)
  seg <- function(v, line){
    d <- with(sol, (v-weight[2]*arange)/weight[1])
    data.frame(t(c(distance=d, angle=arange)), line)
  }
  seg.df <- rbind(seg(1-sol$margin,"margin"),
                  seg(1+sol$margin,"margin"),
                  seg(-1-sol$margin,"margin"),
                  seg(-1+sol$margin,"margin"),
                  seg(1,"decision"),
                  seg(-1,"decision"))
  library(ggplot2)
  comparePlot <- ggplot()+
  geom_point(aes(distance,angle,colour=factor(yi),
                 size=constraint), data=diffs)+
  scale_size_manual(values=c(active=2,inactive=1))+
  geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                   linetype=line),data=seg.df)+
  scale_linetype_manual(values=c(decision="solid",margin="dotted"))
  print(comparePlot)
})

hardCompareQP <- structure(function
### Fit a sparse linear kernel hard margin comparison model to
### linearly separable data. We first normalize the data using
### \code{\link{pairs2svmData}}, resulting in a scaled n x p feature
### difference matrix X, with a new vector of comparisons y in
### c(-1,1). We then define the linear kernel matrix K=XX' and solve
### the dual quadratic program (QP) of SVM: \eqn{\min_{\alpha\in R^n}
### \alpha' K \alpha/2 - y'\alpha} subject to for all i, \eqn{y_i
### \alpha_i \ge 0}. The learned function in the scaled binary SVM
### space is \eqn{f(x) = b + \sum_{i\in sv} \alpha_i k(d_i, x)} where
### sv are the support vectors and the bias b is calculated using the
### average of \eqn{b = y_i - f(d_i)} over all support vectors i. The
### learned ranking function in the original space is \eqn{r(x) =
### \sum_{i\in sv} -\alpha_i/b k(d_i, Sx)} where S is the diagonal
### scaling matrix of the input features. Since we use the linear
### kernel k, we can also write this function as \eqn{r(x) = w'x} with
### the weight vector \eqn{w = -S/b \sum_{i\in sv} \alpha_i d_i}.
(Pairs,
### see \code{\link{check.pairs}}.
 add.to.diag=1e-10,
### This value is added to the diagonal of the kernel matrix, to
### ensure that it is positive definite.
 sv.threshold=1e-3
### Optimal coefficients \eqn{\alpha_i} with absolute value greater
### than this value are considered support vectors.
 ){
  svmData <- pairs2svmData(Pairs)
  sigma <- svmData$scale
  X <- svmData$features
  y <- svmData$labels
  N <- nrow(X)
  P <- ncol(X)
  ## Calculcate the linear kernel and add to diag to make sure it is
  ## numerically positive definite.
  K <- X %*% t(X)
  diag(K) <- diag(K) + add.to.diag
  ## Construct the QP constraints using the quadmod modeling language.
  vars <- make.ids(alpha=N)
  constraints <-
    c(vars$alpha[]*y >= 0,
      list(sum(vars$alpha) == 0))
  ## Fit the linear kernel hard margin SVM using a QP solver.
  sol <- run.quadprog(vars, K, y, constraints)
  sol$sigma <- sigma
  is.sv <- abs(sol$alpha) > sv.threshold
  a.sv <- sol$alpha[is.sv]
  X.sv <- X[is.sv,]
  y.sv <- y[is.sv]
  sol$sv <- list(X=X.sv, a=a.sv, y=y.sv)
  sol$check <- function(X){
    stopifnot(ncol(X)==P)
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))
  }
  weight.svm <- colSums(a.sv * X.sv) 
  bias.values <- y.sv - X.sv %*% weight.svm
  bias <- mean(bias.values)
  sol$margin <- -1/bias
  sol$weight <- -weight.svm/sigma/bias
  sol$rank <- function(X){
    sol$check(X)
    X %*% sol$weight
  }
  sol$predict <- function(Xi, Xip){
    for(X in list(Xi, Xip)){
      sol$check(X)
    }
    D <- Xip - Xi
    rank.diff <- sol$rank(D)
    ifelse(rank.diff < -1, -1L,
           ifelse(rank.diff > 1, 1L, 0L))
  }
  sol
### Comparison model fit. You can do fit$rank(X) to get m numeric
### ranks for the rows of the m x p numeric matrix X. For two feature
### vectors xi and xip, we predict no significant difference if their
### absolute rank difference is less than 1. You can do
### fit$predict(Xi,Xip) to get m predicted comparisons in c(-1,0,1),
### for m by p numeric matrices Xi and Xip. Also, fit$sigma are the
### scales of the input features, fit$sv are the support vectors (in
### the scaled space) and fit$weight is the optimal weight vector (in
### the original space), and if fit$margin is positive than the data
### are separable.
},ex=function(){
  data(separable)
  sol <- hardCompareQP(separable)
  ## check to make sure we have perfect prediction.
  y.hat <- with(separable, sol$predict(Xi, Xip))
  stopifnot(separable$yi == y.hat)
  ## This should also be the same:
  fxdiff <- with(separable, sol$rank(Xip)-sol$rank(Xi))
  y.hat2 <- ifelse(fxdiff < -1, -1L,
                   ifelse(fxdiff > 1, 1L, 0L))
  stopifnot(y.hat == y.hat2)
  ## difference vectors and support vectors to plot.
  point.df <- with(separable, data.frame(Xip-Xi, yi))
  sv.df <- with(sol$sv, data.frame(t(t(X)*sol$sigma)))
  ## calc svm decision boundary and margin.
  mu <- sol$margin
  arange <- range(point.df$angle)
  seg <- function(v, line){
    d <- (v-sol$weight[2]*arange)/sol$weight[1]
    data.frame(t(c(distance=d, angle=arange)), line)
  }
  seg.df <- rbind(seg(1-mu,"margin"),
                  seg(1+mu,"margin"),
                  seg(-1-mu,"margin"),
                  seg(-1+mu,"margin"),
                  seg(1,"decision"),
                  seg(-1,"decision"))
  library(ggplot2)
  svplot <- ggplot()+
  geom_point(aes(distance,angle,colour=factor(yi)), data=point.df,
             size=3)+
  geom_point(aes(distance,angle), data=sv.df,size=1.5)+
  geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                   linetype=line),data=seg.df)+
  scale_linetype_manual(values=c(decision="solid",margin="dotted"))+
  ggtitle(paste("Hard margin linear kernel comparison model",
                "support vectors in black",sep="\n"))
  print(svplot)
})


