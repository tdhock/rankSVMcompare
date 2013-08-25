\name{hardCompareLP}
\alias{hardCompareLP}
\title{hardCompareLP}
\description{Fit a linear hard margin comparison model to linearly separable
data. The LP is max_{\eqn{\mu\in R}, \eqn{w\in R^p}}
mu subject to these constraints: \itemize{ \item for all i such
that \eqn{y_i=0}, \eqn{\mu \le 1-w'(x_i'-x_i)} and \eqn{\mu \le
1+w'(x_i'-x_i)}. \item for all i such that \eqn{y_i\ne 0},
\eqn{\mu \le -1 + w'(x_i'-x_i)y_i}.}}
\usage{hardCompareLP(Pairs)}
\arguments{
  \item{Pairs}{see \code{\link{check.pairs}}.}
}

\value{Comparison model fit. You can do fit$rank(X) to get m numeric
ranks for the rows of the m x p numeric matrix X. For two feature
vectors xi and xip, we predict no significant difference if their
rank difference is less than 1. You can do fit$predict(Xi,Xip) to
get m predicted comparisons in c(-1,0,1), for m by p numeric
matrices Xi and Xip. Also, fit$weight is the optimal vector of p
numeric weights, and if fit$margin is positive then the data are
separable.}

\author{Toby Dylan Hocking}




\examples{
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
}
