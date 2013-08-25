\name{hardCompareQP}
\alias{hardCompareQP}
\title{hardCompareQP}
\description{Fit a sparse linear kernel hard margin comparison model to
linearly separable data. We first normalize the data using
\code{\link{pairs2svmData}}, resulting in a scaled n x p feature
difference matrix X, with a new vector of comparisons y in
c(-1,1). We then define the linear kernel matrix K=DD' and solve
the dual quadratic program (QP) of SVM: \eqn{\min_{\alpha\in R^n}
\alpha' K \alpha/2 - y'\alpha} subject to for all i, \eqn{y_i
\alpha_i \ge 0}. The learned function in the scaled binary SVM
space is \eqn{f(x) = b + \sum_{i\in sv} \alpha_i k(d_i, x)} where
sv are the support vectors and the bias b is calculated using the
average of \eqn{b = y_i - f(d_i)} over all support vectors i. The
learned ranking function in the original space is \eqn{r(x) =
\sum_{i\in sv} -\alpha_i/b k(d_i, Sx)} where S is the diagonal
scaling matrix of the input features. Since we use the linear
kernel k, we can also write this function as \eqn{r(x) = w'x} with
the weight vector \eqn{w = -S/b \sum_{i\in sv} \alpha_i d_i}.}
\usage{hardCompareQP(Pairs, add.to.diag = 1e-10, sv.threshold = 0.001)}
\arguments{
  \item{Pairs}{see \code{\link{check.pairs}}.}
  \item{add.to.diag}{This value is added to the diagonal of the kernel matrix, to
ensure that it is positive definite.}
  \item{sv.threshold}{Optimal coefficients \eqn{\alpha_i} with absolute value greater
than this value are considered support vectors.}
}

\value{Comparison model fit. You can do fit$rank(X) to get m numeric
ranks for the rows of the m x p numeric matrix X. For two feature
vectors xi and xip, we predict no significant difference if their
absolute rank difference is less than 1. You can do
fit$predict(Xi,Xip) to get m predicted comparisons in c(-1,0,1),
for m by p numeric matrices Xi and Xip. Also, fit$sigma are the
scales of the input features, fit$sv are the support vectors (in
the scaled space) and fit$weight is the optimal weight vector (in
the original space), and if fit$margin is positive than the data
are separable.}

\author{Toby Dylan Hocking}




\examples{
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
}
