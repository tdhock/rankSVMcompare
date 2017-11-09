compute_threshold <- function(rankdiff_vector,labls){
  rd<- data.table(diff=rankdiff_vector, label=labls)
  rd[, abs.diff := abs(diff)]
  ord <- rd[order(-abs.diff)]
  ord[, FP := cumsum(label==0)]
  ord[, TP := cumsum((label==-1 & diff <  0) | (label == 1 & diff > 0))]
  ord[, TPR := TP / sum(label != 0)]
  ord[, FPR := FP / sum(label != 0)]
  ord <- rbind(data.table(diff=0,label=0,abs.diff=0,TP=0,TPR=0,FP=0,FPR=0),ord)
  return(rd)
}
