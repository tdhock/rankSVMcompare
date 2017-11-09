library(rankSVMcompare)

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

rankSVM_CV <-
  function(Xi_train,
           Xip_train,
           train_yi,
           parameter.df,
           kern,
           splits_train) {
    
    #Cross validation loop using different kernels
    #The Xi_train,Xip_train and train_yi arguments are ordered by event date,  First element is oldest event, last element is newest
    if(is.character(kern)){
      kern.fun <- get(kern)
    }
    else if(is.function(kern)){
      kern.fun <- kern
    }
    else{
      stop("kern must be either character or function")
    }
    
    if(is.null(parameter.df$C)){
      stop("C must not be null")
    }
    
    store_AUC <- list()
    
    c.col.index <- which(names(parameter.df)=="C")
    
    for (combination.i in 1:nrow(parameter.df)) {
      
      one.combination <- parameter.df[combination.i,]
      one.combination <- data.frame(one.combination)
      
      kern.params <- one.combination[,-c.col.index]
      kern.params <- list(kern.params)
      C.param <- one.combination[,c.col.index]
      svm_kern <- do.call(kern.fun, kern.params)
      
      for (split.i in 1:length(splits_train)) {
        
        rows_train <- 1:(as.integer(nrow(Xi_train) * splits_train[split.i]))
        
        rows_validation <- (as.integer(nrow(Xi_train) * splits_train[split.i]) + 1):as.integer(nrow(Xi_train))
        Xi_train_run <-
          Xi_train[rows_train,]
        Xi_validation_run <-
          Xi_train[rows_validation,]
        
        Xip_train_run <-
          Xi_train[rows_train,]
        Xip_validation_run <-
          Xip_train[rows_validation,]
        train_yi <- data.frame(train_yi)
        train_yi_run <-
          train_yi[rows_train,1]
        validation_yi_run <-
          train_yi[rows_validation,1]
        
        train_pairs <-
          list(
            Xi = Xi_train_run,
            Xip = Xip_train_run,
            yi = train_yi_run
          )
        
        Xi_validation_run[is.na(Xi_validation_run)] <- 0
        Xip_validation_run[is.na(Xip_validation_run)] <- 0
        validation_yi_run[is.na(validation_yi_run)] <- 0
        
        qp <- softCompareQP(train_pairs,svm_kern,C=C.param)
        rankdiff <- rankdiff_qp(qp,Xi_validation_run,Xip_validation_run)
        
        AUC <- calc_AUC(rankdiff,validation_yi_run)

       store_AUC[[paste(combination.i, split.i)]] <-
          data.frame(
            TYPE = paste(one.combination, collapse = "_"),
            CV = splits_train[split.i],
            AUC = AUC
          )
        print("Done1")
      }
    }
    
    store_AUC <- do.call(rbind, store_AUC)
    return(store_AUC)
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

