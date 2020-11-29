processCOVIDData <- function(ntree = 500,nIter = 5, num_folds = 10, num.inner.folds = 7, num.boruta.folds = 5) {
  
  library(pROC)
  library(caret)
  library(gridExtra)
  library(gtable)
  library(ggpubr)
  library(cutpointr)


  setwd('C:/DATI/Policlinico-Humanitas-COVIRT/AARcode/')
  source(file.path('.', 'Utils_preprocess.R'))
  source(file.path('.', 'Utils_postprocess.R'))  
  source(file.path('.', 'boruta_n_fold_select.R'))  
  source(file.path('.','perf.meas.R'))
  source(file.path('.','p.value.analysis.R'))
 
  source(file.path('.','check_MCAR.R'))
  
  
  conf.level = 0.05
  ntree = 12
  num_folds = 10
  nIter.imp = 25
  nIter.class = 5
  n.boot.runs = 100
  
  num.boruta.folds = 10
  dirdata = file.path('data')
  
  
  data = read.csv(file = file.path(dirdata, 'dataEle.csv'),  
                  header = TRUE, sep = ",")
  
  print(paste(1:length(colnames(data)), colnames(data)))
 # LABEL.idx <- readline(prompt="Enter idx of LABEL column: ")
  LABEL.idx <- 6
  
  #start_idx <- as.integer(readline(prompt="Enter starting column data: "))

  start_idx <- 8
  
  LABEL = as.factor(data[,colnames(data)[LABEL.idx]])
  data = data[,colnames(data)[start_idx:ncol(data)]]
  
  
  # check_MCAR pulisce togliende variabili con troppa poca variabilita' ed eventualmente imputando
     data.imp = check_MCAR(data=data,LABEL=LABEL, nIt = nIter.imp)
     df_s = p.value.analysis(data.imp, LABEL = LABEL, n.boot.runs = n.boot.runs)
     write.csv(df_s, file.path(dirdata, 'report_imputed.csv'))
    # 
    # 
     df_s = p.value.analysis(data = data, LABEL = LABEL, n.boot.runs = n.boot.runs)
     write.csv(df_s, file.path(dirdata, 'report_noimpute.csv'))
    # 

    
     
    folds.list = vector("list",length = nIter.class)
    
  
    for ( n.it in 1:nIter.class){
      ###########
      folds.list[[n.it]] <- caret::createFolds(y= LABEL, k = num_folds, list = TRUE, returnTrain = FALSE)
      ##########
    }
  
  
 
  data = data.imp
  num_classes = length(unique(LABEL))
    
  s=1
  
  metrics <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "F1", "Balanced Accuracy" )
  all.metrics = c(metrics, "mcc", "auc")
  results = data.frame(matrix(0, nrow = nIter.class, ncol= length(all.metrics)))
  colnames(results) = all.metrics
  
  for (n.it in 1:nIter.class){
  
    n.it = 1
    ext.folds = folds.list[[n.it]]
    importances = data.frame(matrix(0, nrow = ncol(data), ncol = num_folds))
    rownames(importances) = colnames(data)
    preds = vector(length = nrow(data))
    preds.prob = preds
    for (num_f in 1:num_folds){
      ind.test = ext.folds[[num_f]]
     
      data.train = data[-ind.test ,]
      lab.train = LABEL[-ind.test ]
      cat(sum(lab.train == levels(lab.train)[1])/sum(lab.train ==levels(lab.train)[2]), '\n')
      data.test = data[ind.test ,]
      lab.test = LABEL[ind.test ]

      # Boruta with folds on training set
      votes <- boruta_n_fold_select(data.train, lab.train, num.boruta.folds = num.boruta.folds, select.tentative = TRUE)
      
      idx_del_boruta <- which(votes < 1)
      if ((length(idx_del_boruta)>0) & (length(idx_del_boruta < ncol(data.train)))){
        idx_del = idx_del_boruta
      }else{
        warning("ALL VARIABLES REJECTED BY BURUTA ON ", num.boruta.folds, " FOLDS ");
        print('keep only vars with p-value < 0.01\n')
        df_s = p.value.analysis(data.train, LABEL = lab.train, n.boot.runs = 1)
        # le colonne dei p value sono in colonna 4 per le categoriche e in colonna 4 per le numeriche
        ps = data.frame(df_s[,colnames(df_s)[4]], row.names = rownames(df_s), col.names= 'p-value')
        if (any(ps < 0.01)){
          idx_del = ps >=0.01
        }else{
          print('no p-value  < 0.01: I will keep all the vars: cross your fingers!\n')
          idx_del = NULL
        }
        
      }
      
      cat( num_f, ') delete ', length( idx_del) , ' columns \n')
      cnamesKeep= colnames(data.train)
      cat('the following columns will be used : ', cnamesKeep[-idx_del], '\n')
      if (length(idx_del)>0){
        data.train = data.train[,-idx_del]
        if (length(ind.test)==1){ data.test = data.test[-idx_del] } 
        else{ data.test = data.test[,-idx_del]}
      }
   
      
      clw <- rep(1 / length(lab.train),num_classes)
      names(clw) = levels(as.factor(lab.train))
      spsz <- rep(min(table(as.factor(lab.train))),num_classes)
      names(spsz) = levels(as.factor(lab.train))
      for (nc in levels(as.factor(lab.train))){ clw[nc] = clw[nc]*length(which(as.factor(lab.train) == nc))}
      
      rf_classifier = tuneRF(y = lab.train, x=data.train,  importance=TRUE, classwt = clw,
                             strata = lab.train, sampsize = spsz, ntreeTry=ntree,
                             mtryStart = ncol(data.train)/2, stepFactor=1, improve=0.01,trace=FALSE, plot=FALSE, doBest=TRUE)
      
      
      
      
      preds[ind.test] <- as.integer(predict(rf_classifier,data.test))-1
      preds.prob[ind.test] <-predict(rf_classifier,data.test, type = "prob")[,2]
      
      
      
      imps <-  importance(rf_classifier, type = 1 )
      if (any(imps<0)){ 
        imps[which(imps<0)]=0
      }
      importances[colnames(data.train), num_f] <- imps
    
      
      
    }# end for num_Folds
   
    preds = as.factor(preds)
    
    res.mc = multiclass.confusion.matrix(predictions = preds,ref.labels=LABEL, metrics=metrics)
    results[n.it,] = c(res.mc, "auc" = pROC::auc(LABEL,preds.prob))
    
    
  }# end for nIter
  
  
  
  
  
  
}