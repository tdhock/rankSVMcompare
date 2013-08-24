\name{hardCompareQP}
\alias{hardCompareQP}
\title{hardCompareQP}
\description{Fit a sparse hard margin comparison model to linearly separable
data. The QP is ... TODO.}
\usage{hardCompareQP(Pairs, epsilon = 1e-06)}
\arguments{
  \item{Pairs}{see check.pairs.}
  \item{epsilon}{
}
}



\author{Toby Dylan Hocking}




\examples{
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
}
