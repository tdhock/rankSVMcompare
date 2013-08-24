\name{hardCompareLP}
\alias{hardCompareLP}
\title{hardCompareLP}
\description{Fit a hard margin comparison model to linearly separable data. The
LP is ... TODO!}
\usage{hardCompareLP(Pairs)}
\arguments{
  \item{Pairs}{see check.pairs.}
}



\author{Toby Dylan Hocking}




\examples{
fxdiff <- as.matrix(diffs[,-1]) \%*\% sol$weight
thresh <- function(x)ifelse(x>1,1,ifelse(abs(x)<1,0,-1))
## check to make sure we have perfect prediction.
stopifnot(thresh(fxdiff) == pairs$yi)
margin <- ifelse(pairs$yi==0,{
  1-abs(fxdiff)
},{
  -1 + pairs$yi * fxdiff
})
on.margin <- abs(margin - sol$margin)<1e-6
diffs$constraint <- ifelse(on.margin, "active", "inactive")

## LP Plotting code
slope <- with(sol, -weight[1]/weight[2])
linedf <- function(line, intercept){
  data.frame(slope, intercept, line)
}
line.df <-
  rbind(linedf("decision",c(-1,1)/sol$weight[2]),
        linedf("margin",(c(-1,1,-1,1)*sol$margin+c(1,1,-1,-1))/sol$weight[2]))
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
}
