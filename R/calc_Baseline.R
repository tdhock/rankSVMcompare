calc_Baseline <- function(labls) {
  most_common <- 1
  if (length(labls[labls == -1]) > length(labls[labls == 1])) {
    most_common <- -1
  }
  ROC <- data.frame(FPR=c(0,1),TPR=c(0,length(labls[labls==most_common])/length(labls[labls==-1 | labls == 1])))
  AUC <- trapz(x=ROC$FPR,y=ROC$TPR)
  return(AUC)
}
