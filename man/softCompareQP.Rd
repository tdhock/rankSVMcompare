\name{softCompareQP}
\alias{softCompareQP}
\title{softCompareQP}
\description{Fit a soft-margin comparison model by using \code{\link{ksvm}}
(libsvm) to solve the dual SVM problem. TODO: explain optimization
problem.}
\usage{softCompareQP(Pairs, ...)}
\arguments{
  \item{Pairs}{see \code{\link{check.pairs}}.}
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
pmodel <- p+
  geom_contour(aes(distance, angle, z=f), size=1.5,
               data=X.grid, colour="grey")+
  geom_dl(aes(distance, angle, z=f, label=..level..), colour="grey",
          data=X.grid, method="bottom.pieces", stat="contour")+
  geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                   linetype=line),data=seg.df)+
  geom_point(aes(distance, angle, shape=sv.type), data=sv.df)+
  scale_linetype_manual(values=c(margin="dashed",decision="solid"))+
  scale_shape_manual(values=c(margin=13,error=3))
print(pmodel)
}
