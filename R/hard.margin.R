hardCompareLP <- structure(function
### Fit a hard margin comparison model to linearly separable data. The
### LP is ... TODO!
(Pairs
### see check.pairs.
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
  result$rank <- function(X){
    stopifnot(ncol(X)==P)
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))
    X %*% result$weight
  }
  result$predict <- function(Xi, Xip){
    for(X in list(Xi, Xip)){
      stopifnot(ncol(X)==P)
      stopifnot(is.matrix(X))
      stopifnot(is.numeric(X))
    }
    D <- Xip-Xi
    rank.diff <- result$rank(D)
    ifelse(rank.diff < -1, -1L,
           ifelse(rank.diff > 1, 1L, 0L))
  }
  result
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
### Fit a sparse hard margin comparison model to linearly separable
### data. The QP is ... TODO.
(Pairs,
### see check.pairs.
 epsilon=1e-6
 ){
  check.pairs(Pairs)
  yi <- Pairs$yi
  is.zero <- yi == 0
  Di <- with(Pairs, scale(Xip - Xi)) ## NEED TO SCALE HERE!!
  Di[!is.zero,] <- Di[!is.zero,] * yi[!is.zero]
  yi[!is.zero] <- 1
  Di.other <- Di
  Di.other[is.zero,] <- -Di.other[is.zero,]
  Di.both <- rbind(Di,Di.other[is.zero,])
  yi.both <- c(yi, rep(0, sum(is.zero)))
  svm.y <- ifelse(pairs$yi==0, -1, 1)
  X <- as.matrix(pairs[,c("distance","angle")])
  K <- X %*% t(X)
  N <- nrow(X)
  vars <- make.ids(alpha=N)
  constraints <-
    c(vars$alpha[]*svm.y >= 0,
      list(sum(vars$alpha) == 0))
  diag(K) <- diag(K)+epsilon
  sol <- run.quadprog(vars, K, svm.y, constraints)
  is.sv <- abs(sol$alpha)>epsilon ## THIS IS BIG!!
  a.sv <- sol$alpha[is.sv]
  X.sv <- X[is.sv,]
  y.sv <- svm.y[is.sv]
  sv.evals <- (a.sv * X.sv) %*% t(X.sv)
  bias.values <- y.sv - colSums(sv.evals)
  bias <- mean(bias.values)
  sol
},ex=function(){
  ## calc svm decision boundary and margin.
  weight <- colSums(X.sv * a.sv / -bias)
  mu <- -1/bias
  arange <- range(diffs$angle)
  seg <- function(v, line){
    d <- (v-weight[2]*arange)/weight[1]
    data.frame(t(c(distance=d, angle=arange)), line)
  }
  seg.df <- rbind(seg(1-mu,"margin"),
                  seg(1+mu,"margin"),
                  seg(-1-mu,"margin"),
                  seg(-1+mu,"margin"),
                  seg(1,"decision"),
                  seg(-1,"decision"))
  model.segs <- rbind(model.segs, {
    data.frame(seg.df, set.name, model="svm")
  })
  model.sv <- rbind(model.sv, data.frame(X.sv, set.name))

  mplot <- ggplot()+
  geom_point(aes(distance,angle,colour=factor(yi),
                 size=constraint), data=model.points)+
  scale_colour_manual(values=yi.colors)+
  scale_size_manual(values=c(active=2,inactive=1))+
  geom_abline(aes(slope=slope,intercept=intercept,linetype=line),
              data=model.lines)+
  geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                   linetype=line),data=model.segs, colour="red")+
  scale_linetype_manual(values=c(decision="solid",margin="dotted"))+
  facet_grid(.~set.name)+
  theme_bw()+
  geom_point(aes(distance, angle), data=model.sv, size=1)+
  theme(panel.margin=unit(0,"cm"))

print(mplot)
})


