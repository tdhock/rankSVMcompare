rankSVM_CV <-
  function(Xi_train,
           Xip_train,
           train_yi,
           parameter.df,
           kern,
           splits_train) {
    
    #Cross validation loop using different kernels
    #The Xi_train,Xip_train and train_yi arguments are ordered by event date,  First element is oldest event, last element is newest
    
    kern.fun <- if(is.character(kern)){
      get(kern)
    }else if(is.function(kern)){
      kern.fun
    }else{
      stop("kern must be either character or function")
    }
    
    if(is.null(parameter.df$C)){
      stop("C must not be null")
    }
  
    store_AUC <- list()
    
    
    c.col.index <- which(names(parameter.df)=="C")
    
    for (combination.i in 1:nrow(parameter.df)) {
      one.combination <- parameter.df[combination.i, ]
      kern.params <- one.combination[, -c.col.index]
      C.param <- one.combination[["C"]]
      svm_kern <- do.call(kern.fun, kern.params)
      for (split.i in 1:length(splits_train)) {
        
        rows_train <- 1:(as.integer(nrow(Xi_train) * splits_train[split.i]))
        rows_validation <- -rows_train
        Xi_train_run <-
          Xi_train[rows_train,]
        Xi_validation_run <-
          Xi_train[rows_validation,]
        
        Xip_train_run <-
          Xi_train[rows_train,]
        Xip_validation_run <-
          Xip_train[rows_validation,]
        
        train_yi_run <-
          train_yi[rows_train,]
        validation_yi_run <-
          train_yi[rows_validation,]
        
        train_pairs <-
          list(
            Xi = Xi_train_run,
            Xip = Xip_train_run,
            yi = train_yi_run
          )
        
        qp <- softCompareQP(train_pairs,svm_kern)
        
        Xi_rank <- rankdiff(qp, Xi_validation_run, Xip_validation_run)
        AUC <- calc_AUC(Xi_rank, validation_yi_run)
        
        store_AUC[[paste(combination.i, split.i)]] <-
          data.frame(
            TYPE = paste(one.combination, collapse = "_"),
            CV = splits_train[split.i],
            AUC = AUC
          )
      }
    }
    
    store_AUC <- do.call(rbind, store_AUC)
    return(store_AUC)
  }
