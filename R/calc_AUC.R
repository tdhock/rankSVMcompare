require(caTools)

calc_AUC <-
  function(qp,
           test_matrix_Xi,
           test_matrix_Xip,
           test_yi) {
    rankdiff <- function(qp, matrix_Xi, matrix_Xip, t) {
      Xirank <- qp$rank.scaled(X = matrix_Xi)
      Xiprank <- qp$rank.scaled(X = matrix_Xip)
      Xirank <- cbind(Xirank, Xiprank)
      
      Xirank <- data.frame(Xirank)
      
      for (i in 1:nrow(Xirank)) {
        diff <- Xirank[i, 2] - Xirank[i, 1]
        if (diff < -1 * t) {
          Xirank[i, 3] <- -1
        }
        if (diff > t) {
          Xirank[i, 3] <- 1
        }
        if (abs(diff) <= t) {
          Xirank[i, 3] <- 0
        }
      }
      return(Xirank)
    }
    
    out_AUC <- data.frame(tau = 0,
                          TP = 0,
                          FP = 0,
                          i = 0)
    
    rankdiff_out <-
      rankdiff(qp, test_matrix_Xi, test_matrix_Xip,t = 1)
    rankdiff_out$difference <- abs(rankdiff_out$X2 - rankdiff_out$X1)
    
    tau_array <- c(abs(rankdiff_out$difference))
    
    out_AUC_SVM <- data.frame(tau = 0, TP = 0, FP = 0)
    
    for (i in tau_array) {
      rankdiff_out <-
        rankdiff(qp, test_matrix_Xi, test_matrix_Xip, i)
      colnames(rankdiff_out)[3] <- "diff"
      rankdiff_out <- cbind(rankdiff_out, test_yi)
      colnames(rankdiff_out)[4] <- "test_yi"
      rankdiff_out$TP <- 0
      rankdiff_out$FP <- 0
      
      #TP
      rankdiff_out[rankdiff_out$diff == rankdiff_out$test_yi &
                     rankdiff_out$test_yi == -1, 5] <- 1
      rankdiff_out[rankdiff_out$diff == rankdiff_out$diff &
                     rankdiff_out$test_yi == 1, 5] <- 1
      
      rankdiff_out[rankdiff_out$diff == -1 &
                     rankdiff_out$test_yi == 0 , 6] <- 1
      rankdiff_out[rankdiff_out$diff == 1 &
                     rankdiff_out$test_yi == 0 , 6] <- 1
      
      TP = nrow(rankdiff_out[rankdiff_out$TP == 1,]) / (nrow(rankdiff_out[rankdiff_out$test_yi == 1,]) + nrow(rankdiff_out[rankdiff_out$test_yi == -1,]))
      FP = nrow(rankdiff_out[rankdiff_out$FP == 1,]) / (nrow(rankdiff_out[rankdiff_out$test_yi == 0,]))
      out_AUC_SVM <-
        rbind(out_AUC_SVM, data.frame(tau = i, TP = TP, FP = FP))
    }
    
    out_AUC_SVM <- out_AUC_SVM[-1, ]
    out_AUC_SVM$tau <- NULL
    yi <- data.frame(yi)
    colnames(yi) <- "matches.RESULT"
    AUC_baseline <-
      data.frame(TP = c(length(yi[yi$matches.RESULT == -1, ]) / (length(yi[yi$matches.RESULT ==
                                                                             -1, ]) + length(yi[yi$matches.RESULT == 1, ])), 0), FP = c(1, 0 / length(yi[yi$matches.RESULT ==
                                                                                                                                                           0, ])))
    AUC_list <- data.frame(TYPE = 0, AUC = 0)
    out_AUC_SVM <- out_AUC_SVM[order(out_AUC_SVM$TP), ]
    AUC <-
      data.frame(TYPE = "TEST",
                 AUC = trapz(x = out_AUC_SVM$TP, y = out_AUC_SVM$FP))
    AUC_baseline <- AUC_baseline[order(AUC_baseline$TP), ]
    Baseline_AUC <-
      data.frame(TYPE = "BASELINE",
                 AUC = trapz(x = AUC_baseline$TP, y = AUC_baseline$FP))
    AUC_list <- rbind(AUC_list, AUC)
    AUC_list <- rbind(AUC_list, Baseline_AUC)
    AUC_list <- AUC_list[-1, ]
    
    return_list <-
      list("AUC_list" = AUC_list,
           "Baseline_ROC" = AUC_baseline,
           "ROC" = AUC)
    
    return(return_list)
    
  }
  