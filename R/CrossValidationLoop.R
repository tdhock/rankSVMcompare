rankSVM_CV <-
  function(Xi_train,
           Xip_train,
           train_yi,
           parameter.df,
           kern) {
    
    #Cross validation loop using different kernels
    #The Xi_train,Xip_train and train_yi arguments are ordered by event date,  First element is oldest event, last element is newest
    
    if (kern == "vanilladot") {
      kern.fun <- kernlab::vanilladot()
    } else if (kern == "rbfdot") {
      kern.fun <- kernlab::rbfdot()
    } else if (kern == "rbfdot") {
      kern.fun <- kernlab::polydot()
    } else{
      return(0)
    }
    
    if(is.null(parameter.df$C)){
      return(0)
    }
  
    store_AUC <- list()
    splits_train <- c(1 / 2, 2 / 3, 3 / 4, 4 / 5, 5 / 6)
    splits_test <- c(1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 6)
    
    for (combination.i in 1:nrow(parameter.df)) {
      for (j in 1:length(splits_test)) {
        one.combination <- parameter.df[combination.i, ]
        Xi_train_run <-
          Xi_train[1:(as.integer(nrow(Xi_train) * splits_train[j])),]
        Xi_test_run <-
          Xi_train[(as.integer(nrow(Xi_train) * splits_train[j] + 1)):(nrow(Xi_train)),]
        
        Xip_train_run <-
          Xi_train[1:(as.integer(nrow(Xi_train) * splits_train[j])),]
        Xip_test_run <-
          Xip_train[(as.integer(nrow(Xip_train) * splits_train[j] + 1)):(nrow(Xip_train)),]
        
        train_yi_run <-
          train_yi[1:(as.integer(nrow(train_yi) * splits_train[j])),]
        test_yi_run <-
          train_yi[(as.integer(nrow(train_yi) * splits_train[j] + 1)):(nrow(train_yi)),]
        
        train_pairs <-
          list(
            Xi = Xi_train_run,
            Xip = Xip_train_run,
            yi = train_yi_run
          )
        test_pairs <-
          list(
            Xi = Xi_test_run,
            Xip = Xip_test_run,
            yi = test_yi_run
          )
        
        qp <- do.call(softCompareQP,args=c(Pairs=train_pairs,one.combination,kernel=kern))
        
        Xi_rank <- rankdiff(qp, Xi_test_run, Xip_test_run)
        AUC <- calc_AUC(Xi_rank, test_yi_run)
        
        store_AUC[[paste(i, j)]] <-
          data.frame(
            TYPE = paste(one.combination, collapse = "_"),
            CV = splits_train[j],
            AUC = AUC
          )
      }
    }
    
    store_AUC <- do.call(rbind, store_AUC)
    return(store_AUC)
  }
