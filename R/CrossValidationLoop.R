rankSVM_CV <-
  function(Xi_train,
           Xip_train,
           test,
           train_yi,
           test_yi,
           C_vector,
           kern,
           sigma_vector) {
    
    #Cross validation loop using different kernels
    
    if (kern == "vanilladot") {
      store_AUC <- list()
      splits_train <- c(1 / 2, 2 / 3, 3 / 4, 4 / 5, 5 / 6)
      splits_test <- c(1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6)
      for (i in 1:length(C_vector)) {
        for (j in 1:length(splits_test)) {
          Xi_train_run <- Xi_train[1:(nrow(Xi_train * splits_train[j])), ]
          Xi_test_run <-
            Xi_train[(nrow(Xi_train * splits_train[j]) + 1):(nrow(Xi_train)), ]
          
          Xip_train_run <-
            train[1:(nrow(train * splits_train[j])), ]
          Xip_test_run <-
            train[(nrow(train * splits_train[j]) + 1):(nrow(train)), ]
          
          train_yi_run <-
            train_yi[1:(nrow(train_yi * splits_train[j])), ]
          test_yi_run <-
            train_yi[(nrow(train_yi * splits_train[j]) + 1):(nrow(train_yi)), ]
          
          train_pairs <-
            list(
              Xi = Xi_train_run,
              Xip = Xip = Xip_train_run,
              yi = train_yi_run
            )
          test_pairs <-
            list(
              Xi = Xi_test_run,
              Xip = Xip = Xip_test_run,
              yi = test_yi_run
            )
          
          qp <-
            softCompareQP(pairs, C = C_vector[[i]], kernel = kern)
          
          Xi_rank <- rankdiff(qp, Xi_test_run, Xip_test_run)
          AUC <- calc_AUC(Xi_rank, test_yi_run)
          
          store_AUC[[i]] <-
            data.frame(TYPE = C_vector[[i]],
                      CV=splits_train[j],
                      AUC = AUC)
        }
      }
    }
    if (kern == "rbfdot" | kern="polydot") {
      store_AUC <- list()
      splits_train <- c(1 / 2, 2 / 3, 3 / 4, 4 / 5, 5 / 6)
      splits_test <- c(1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6)
      for (k in 1:length(sigma_vector)) {
        for (i in 1:length(C_vector)) {
          for (j in 1:length(splits_test)) {
            Xi_train_run <- Xi_train[1:(nrow(Xi_train * splits_train[j])), ]
            Xi_test_run <-
              Xi_train[(nrow(Xi_train * splits_train[j]) + 1):(nrow(Xi_train)), ]
            
            Xip_train_run <-
              train[1:(nrow(train * splits_train[j])), ]
            Xip_test_run <-
              train[(nrow(train * splits_train[j]) + 1):(nrow(train)), ]
            
            train_yi_run <-
              train_yi[1:(nrow(train_yi * splits_train[j])), ]
            test_yi_run <-
              train_yi[(nrow(train_yi * splits_train[j]) + 1):(nrow(train_yi)), ]
            
            train_pairs <-
              list(
                Xi = Xi_train_run,
                Xip = Xip = Xip_train_run,
                yi = train_yi_run
              )
            test_pairs <-
              list(
                Xi = Xi_test_run,
                Xip = Xip = Xip_test_run,
                yi = test_yi_run
              )
            
            qp <-
              softCompareQP(pairs,
                            C = C_vector[[i]],
                            kernel = kern,
                            sigma = sigma_vector[[k]])
            
            Xi_rank <- rankdiff(qp, Xi_test_run, Xip_test_run)
            AUC <- calc_AUC(Xi_rank, test_yi_run)
            
            store_AUC[[i]] <-
              data.frame(TYPE = paste(sigma_vector[[k]],C_vector[[i]], sep = "_"),
                         CV=splits_train[j],
                         AUC = AUC)
          }
        }
      }
    }
    else{
      return(-1)
    }
    store_AUC <- do.call(rbind, store_AUC)
    return(store_AUC)
  }
