compute_FP_threshold <- function(rankdiff_vector,labls){
  rd<- data.table(diff=rankdiff_vector, label=labls)
  rd[, abs.diff := abs(diff)]
  rd[order(abs.diff)]
  rd[, FP := cumsum(label==0)]
  rd[, FPR := FP / sum(label != 0)]
  return(rd)
}

compute_TP_threshold <- function(rankdiff_vector,labls){
  rd<- data.table(diff=rankdiff_vector, label=labls)
  rd[, abs.diff := abs(diff)]
  rd[order(abs.diff)]
  rd[, TP := cumsum((label==-1 & diff <  0) | (label == 1 & diff > 0))]
  rd[, TPR := TP / sum(label != 0)]
  rd <- rbind(data.table(diff=0,label=0,abs.diff=0,TP=0,TPR=0),rd)
  return(rd)
}
