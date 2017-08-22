require(caTools)
require(PlayerRatings)

compute_AUC_andELO <- function(train, test) {
  
## Require chess match data format for train and test df colnames(Month, Whote, Black, Score)  
  
  temp <- tempfile()
  download.file(paste(cm, "initial_ratings.csv", sep = ""), temp)
  cSt <- read.csv(temp)
  cSt <-
    data.frame(
      Player = cSt$Player,
      Rating = 1200,
      Deviation = 200,
      Games = 0,
      Win = 0,
      Draw = 0,
      Loss = 0,
      Lag = 0
    )
  
  robjf2 <- fide(train, cSt)
  re <- robjf2$ratings
  cSt <- re[, c()]
  cSt <- re[, c(1, 2, 3, 4, 5, 6, 7)]
  cSt$Deviation <- 200
  cSt <- cSt[, c(1, 2, 8, 3, 4, 5, 6, 7)]
  
  pred_array <- list()  
  array_counter <- 1
  
  for (i in 1:nrow(test)) {
    prf <- predict(robjf2, test[i,], tng = 0)
    robjf2 <- fide(test[i,], cSt)
    re <- robjf2$ratings
    cSt <- re[, c(1, 2, 3, 4, 5, 6, 7)]
    cSt$Deviation <- 200
    cSt <- cSt[, c(1, 2, 8, 3, 4, 5, 6, 7)]
    
    pred_array[[array_counter]] <- data.frame(pred = prf)
    array_counter <- array_counter + 1
  }
  pred_array <- do.call(rbind, pred_array)    

  ##AUC ROC
  
  test <- cbind(test, pred_array)
  
  test$Month <- NULL
  test$White <- NULL
  test$Black <- NULL
  
  tau_array <- c(abs(0.5 - test$pred))

  out_AUC <- list()
  array_counter <- 1
  for (i in tau_array) {
    test$pred_bounded <- 0
    range_top <- 0.5 + i
    range_bottom <- 0.5 - i
    for (j in 1:nrow(test)) {
      if (test[j, 2] >= range_bottom & test[j, 2] <= range_top) {
        test[j, 3] <- 0.5
      }
      else {
        if (test[j, 2] < range_bottom) {
          test[j, 3] <- 0
        }
        if (test[j, 2] > range_top) {
          test[j, 3] <- 1
        }
      }
    }
    
    test$TP <- 0
    test$FP <- 0
    
    test[test$Score == 0 & test$pred_bounded == 0, 4] <- 1
    test[test$Score == 1 & test$pred_bounded == 1, 4] <- 1
    
    test[test$Score == 0.5 & test$pred_bounded == 0, 5] <- 1
    test[test$Score == 0.5 & test$pred_bounded == 1, 5] <- 1
    
    TP <-
      nrow(test[test$TP == 1,]) / (nrow(test[test$Score == 1,]) + nrow(test[test$Score == 0,]))
    FP <-
      nrow(test[test$FP == 1,]) / nrow(test[test$Score == 0.5,])
    
    out_AUC[[array_counter]] <- data.frame(tau = i, TP = TP, FP = FP)
    array_counter <- array_counter + 1
  }
  
  out_AUC <- do.call(rbind, out_AUC) 
  out_AUC$tau <- NULL
    
  out_AUC <- out_AUC[order(out_AUC$FP,out_AUC$TP),]
  AUC <- data.frame(TYPE="ELO",AUC=caTools::trapz(out_AUC$TP,out_AUC$FP))
  
  AUC_calc <- 0
  for (i in 1:(nrow(out_AUC)-1)){
    AUC_calc <- AUC_calc + (out_AUC$FP[i+1] - out_AUC$FP[i])*out_AUC$TP[i]
  }
  
  return_list <- list("AUC"=AUC,"ROC"=out_AUC)
  return(return_list)
  
}
