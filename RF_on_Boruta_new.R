RF_on_Boruta_new <- function(data, LABEL, ntree = 100, num_folds =5, num.boruta.folds = 3){
  # legge i dati nel file data_filename, i tipi nel file type_filename, i folds nella matrice folds
  # data_types contiene LABELS come tipo delle colonne da usare come labels
  library(pracma)
  cat('Feature selection with inner ', num_folds, '-fold RF with inner ',  num.boruta.folds , '-fold Boruta', '\n')
 
  preds_RF_multi <- as.factor(LABEL)
  data_matrix <- data.matrix(data) 
  
  
  importances = data.frame(matrix(0, nrow = ncol(data), ncol = num_folds))
  rownames(importances) = colnames(data)
  
  num_classes = length(levels(LABEL))
  
  inner.folds <- createFolds(y= LABEL, k = num_folds, list = TRUE, returnTrain = FALSE)
  invotes = rep(0,times = ncol(data))
  
  for (num_f in 1:num_folds){
    cat('internal fold', num_f, '\n')
    ind.test = inner.folds[[num_f]]     
    # per la matrice di importanza attuale uso le colonne dei dati ch ho ora!
   
    lab.train <- LABEL[-ind.test]
    
    cat("In RF Boruta_new LABEL is categorical? ", is.factor(lab.train), '\n')
    data.train <- data_matrix[-ind.test, ]
    lab.test <- LABEL[ind.test]
    data.test <- data_matrix[ind.test,]
    
    
    cat('.', num_f, '\n')
    
    votes <- boruta_n_fold_select(data.train, lab.train, num_folds = num.boruta.folds, select.tentative = TRUE)
    
    idx_del_boruta <- which(votes < pracma::ceil(num.boruta.folds/2))
    if (length(idx_del_boruta)>0){
      data.train = data.train[,-idx_del_boruta]
      data.test = data.test[,-idx_del_boruta]
    }
    if (ncol(data.train)>0){
    
      clw <- rep(1 / length(lab.train),num_classes)
      names(clw) = levels(as.factor(lab.train))
      spsz <- rep(min(table(as.factor(lab.train))),num_classes)
      names(spsz) = levels(as.factor(lab.train))
      for (nc in levels(as.factor(lab.train))){ clw[nc] = clw[nc]*length(which(as.factor(lab.train) == nc))}
       
      rf_classifier = tuneRF(y = lab.train, x=data.train,  importance=TRUE, classwt = clw,
             strata = lab.train, sampsize = spsz, ntreeTry=ntree,
             mtryStart = ncol(data.train)/2, stepFactor=1, improve=0.01,trace=FALSE, plot=FALSE, doBest=TRUE)
      
      
      
      
      preds_RF_multi[ind.test] <- predict(rf_classifier,data.test)
      imps <-  importance(rf_classifier, type = 1 )
      if (any(imps<0)){ 
        imps[which(imps<0)]=0
      }
      importances[colnames(data.train), num_f] <- imps
    }else{
      warning("ALL VARIABLES REJECTED BY BURUTA ON ", num.boruta.folds, " FOLDS");
      
    }
   
  }
 
  result <- list("importances" = importances, "predictions" = preds_RF_multi, 
                 "LABELS" = LABEL )
  return(result)
}