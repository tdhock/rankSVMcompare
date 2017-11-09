library(rankSVMcompare)
library(PlayerRatings)

trapz <- function (x, y)
{
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1])) / 2)
}

rankdiff_qp <- function(qp, matrix_Xi, matrix_Xip, t) {
  Xirank <- qp$rank(X = matrix_Xi)
  Xiprank <- qp$rank(X = matrix_Xip)
  Xirank <- cbind(Xirank, Xiprank)
  
  Xirank <- data.frame(Xirank)
  Xirank$diff <- Xirank[, 2] - Xirank[, 1]
  Xirank$label <- NA
  
  for (i in 1:nrow(Xirank)) {
    diff <- Xirank[i, 2] - Xirank[i, 1]
    if (diff < (-1 * t)) {
      Xirank[i, 4] <- -1
    }
    if (diff > t) {
      Xirank[i, 4] <- 1
    }
    if (abs(diff) <= t) {
      Xirank[i, 4] <- 0
    }
  }
  return(Xirank)
}

rankdiff <- function(rank, rank_p, t) {

  Xirank <- cbind(rank, rank_p)
  
  Xirank <- data.frame(Xirank)
  Xirank$diff <- Xirank[, 2] - Xirank[, 1]
  Xirank$label <- NA
  
  for (i in 1:nrow(Xirank)) {
    diff <- Xirank[i, 2] - Xirank[i, 1]
    if (diff < (-1 * t)) {
      Xirank[i, 4] <- -1
    }
    if (diff > t) {
      Xirank[i, 4] <- 1
    }
    if (abs(diff) <= t) {
      Xirank[i, 4] <- 0
    }
  }
  return(Xirank)
}

compute_FP_threshold <- function(rankdiff_vector, labls) {
  rd <- data.table(diff = rankdiff_vector, label = labls)
  rd[, abs.diff := abs(diff)]
  ord <- rd[order(-abs.diff)]
  ord[, FP := cumsum(label == 0)]
  ord[, FPR := FP / sum(label != 0)]
  return(rd)
}

compute_threshold <- function(rankdiff_vector, labls) {
  rd <- data.table(diff = rankdiff_vector, label = labls)
  rd[, abs.diff := abs(diff)]
  ord <- rd[order(-abs.diff)]
  ord[, FP := cumsum(label == 0)]
  ord[, TP := cumsum((label == -1 &
                        diff <  0) | (label == 1 & diff > 0))]
  ord[, FPR := FP / sum(label != 0)]
  ord[, TPR := TP / sum(label != 0)]
  ord <-
    rbind(data.table(
      diff = 0,
      label = 0,
      abs.diff = 0,
      FP = 0,
      TP = 0,
      FPR = 0,
      TPR = 0
    ),
    ord)
  return(ord)
}

calc_AUC <- function(rankdiff_vector, labls) {
  ROC <- compute_threshold(rankdiff_vector, labls)
  AUC <- trapz(x = ROC$FPR, y = ROC$TPR)
  return(AUC)
}

get_ROC <- function(rankdiff_vector, labls) {
  ROC <- compute_threshold(rankdiff_vector, labls)
  return(ROC)
}

calc_Baseline <- function(labls) {
  most_common <- 1
  if (length(labls[labls == -1]) > length(labls[labls == 1])) {
    most_common <- -1
  }
  ROC <- data.frame(FPR=c(0,1),TPR=c(0,length(labls[labls==most_common])/length(labls[labls==-1 | labls == 1])))
  AUC <- trapz(x=ROC$FPR,y=ROC$TPR)
  return(AUC)
}


## TESTS calcAUC
test_rankdiff_vector <- c(0.1, 0.2, -0.3, 0.4)
test_labels <- c(1, 0, -1, 0)
AUC_test <- calc_AUC(test_rankdiff_vector, test_labels)

if (AUC_test != 0.25) {
  stop()
}

test_rankdiff_vector <- c(0,1)
test_labels <- c(0,1)
AUC_test <- calc_AUC(test_rankdiff_vector, test_labels)
if (AUC_test != 1) {
  stop()
}

# TESTS Baseline
test_labels <- c(1,0,-1,0,1,1)
AUC_test <- calc_Baseline(test_labels)
if(AUC_test != 0.375){
  stop()
}
