trapz <- function (x, y)
{
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx -1]))/2)
}

calc_AUC <- function(rankdiff_vector,labls){
  TPR <- compute_TP_threshold(rankdiff_vector,labls)
  FPR <- compute_FP_threshold(rankdiff_vector,labls)
  AUC <- trapz(x=FPR$FPR,y=TPR$TPR)
  return(AUC)
}
