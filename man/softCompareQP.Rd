\name{softCompareQP}
\alias{softCompareQP}
\title{softCompareQP}
\description{Fit a soft-margin comparison model by using \code{\link{ksvm}}
(libsvm) to solve the dual SVM problem. TODO: explain optimization
problem.}
\usage{softCompareQP(Pairs, kernel = rbfdot(sigma = 1), ...)}
\arguments{
  \item{Pairs}{see \code{\link{check.pairs}}.}
  \item{kernel}{Kernel function, see \code{\link{dots}}.}
  \item{\dots}{Passed to \code{\link{ksvm}}.}
}



\author{Toby Dylan Hocking}




\examples{
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
  outputs <- list(svm.f=fit$svm.f(sc.grid),
                  rank=fit$rank.scaled(sc.grid))
  for(fun.name in names(outputs)){
    scaledGrid <- rbind(scaledGrid,{
      data.frame(cost, sc.grid, value=outputs[[fun.name]], fun.name)
    })
  }
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
      sprintf("C = \%s",as.character(val))
    }else as.character(val)
  })+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  ggtitle("Support vector comparison model (grey level lines)")
print(diffPlot)
## Fit soft-margin non-linear comparison models for some cost and
## kernel width parameters.
grid.df <- data.frame()
for(cost in c(1)){
  for(sigma in c(1/2, 2)){
    fit <- softCompareQP(not, kernel=rbfdot(sigma), C=cost)
    f <- fit$rank(X.grid)
    grid.df <- rbind(grid.df, data.frame(X.grid, f, cost, sigma))
  }
}
library(directlabels)
library(grid)
nonlinear <- arrowPlot+
  geom_contour(aes(distance, angle, z=f), size=1.5,
               data=grid.df, colour="grey")+
  geom_dl(aes(distance, angle, z=f, label=..level..), colour="grey",
          data=grid.df, method="bottom.pieces", stat="contour")+
  facet_grid(sigma~yi,labeller=function(var, val){
    sprintf("\%s = \%s",var,as.character(val))
  })+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  ggtitle("Support vector comparison model (grey level curves)")
print(nonlinear)
}
