main_clinical_data_analysis <- function(ntree = 500,nIter = 5, num_folds = 10, num.inner.folds = 7, num.boruta.folds = 5) {
  
  library(pROC)
  library(caret)
  library(gridExtra)
  library(gtable)
  
  
  conf.level = 0.05
  ntree = 12
  num_folds = 10
  nIter = 5
  num.inner.folds = 7 
  num.boruta.folds = 5
  dirdata = file.path('C:', 'DATI', 'Policlinico-Humanitas-COVIRT', 'Caterina-AI-Covid', 'clinical')
  
  setwd(dirdata)
  
  data = read.csv(file = file.path(dirdata, 
                                   paste('data_ele_noTOT.csv', sep="",collapse = NULL)),  
                  header = TRUE, sep = ",")
  
  LABEL = as.factor(data$LABEL)
  NV = (ncol(data)-1)/2
  cat(length(unique(LABEL)))
  
 
  metrics <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "F1", "Balanced Accuracy" ) 
  folds.list = vector("list",length = nIter)
  for ( n.it in 1:nIter){
    ###########
    folds.list[[n.it]] <- caret::createFolds(y= LABEL, k = num_folds, list = TRUE, returnTrain = FALSE)
    ##########
  }
  
  
 # wilcox.test(auc ~ risk.model, data = results_missForest_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
                                     
  
  meanTabArea = as.data.frame(matrix (0, ncol= NV, nrow= length(unique(LABEL))))
  colnames(meanTabArea) = colnames(data)[1:NV]
  rownames(meanTabArea) = unique(LABEL)
  stderrTabArea = as.data.frame(matrix (0, ncol= NV, nrow= length(unique(LABEL))))
  rownames( stderrTabArea) = unique(LABEL)
  colnames( stderrTabArea) = colnames(data)[1:NV]
  
  meanPercArea =as.data.frame(matrix (0, ncol= NV, nrow= length(unique(LABEL))))
  rownames(meanPercArea) = unique(LABEL)
  colnames(meanPercArea) = colnames(data)[1:NV]
  
  stderrPercArea =as.data.frame(matrix (0, ncol= NV, nrow= length(unique(LABEL))))
  rownames(stderrPercArea) = unique(LABEL)
  colnames(stderrPercArea) = colnames(data)[1:NV]
  cnames=colnames(data)
  cnames= cnames[-grep("LABEL",cnames)]
  wilcoxArea = data.frame(matrix(nrow=3,ncol =NV*2))
  colnames(wilcoxArea) = rbind(paste(cnames[1:NV], "greater"),paste(cnames[1:6], "lower"))
  rownames(wilcoxArea) = list("discharged vs hospitalized", "disharged vs ICU", "hospitalized vs ICU")
  for (nC in 1:NV){
    trueNC = grep(cnames[nC], colnames(data)) 
    cat('column  = ',  colnames(data)[trueNC], 'idx = ', trueNC, '\n')
   
    for (numL in 1:length(unique(LABEL))){
      nL= unique(LABEL)[[numL]]
      meanTabArea[unique(LABEL)==nL,nC] = mean(data[LABEL == nL,trueNC])
      stderrTabArea[unique(LABEL)==nL,nC] = sd(data[LABEL == nL,trueNC])/sqrt(sum(LABEL==nL))
      cat( colnames(data)[trueNC], ' (', nL , ') = ', meanTabArea[unique(LABEL)==nL,nC] , ' +- ', stderrTabArea[unique(LABEL)==nL,nC], '\n' ) 
      
    }
    
    cat('********************wilcoxon test \n')
   
    wilcoxArea[ "disharged vs hospitalized", nC*2-1: nC*2]=      
        cbind('greater' = wilcox.test(x= data[LABEL == 0, trueNC], y = data[LABEL == 1, trueNC],  
                                      paired = FALSE, alternative = "greater")[["p.value"]],
              'lower'=wilcox.test(x= data[LABEL == 0, trueNC], y = data[LABEL ==1, trueNC],  
                                  paired = FALSE, alternative = "greater")[["p.value"]])
    
    
    
    
    trueNC = trueNC+NV
    cat('column  = ',  colnames(data)[trueNC], '\n')
    
    
    
    for (nL in unique(LABEL)){
           # for (nL_other in unique(LABEL)){
    #    if (nL_other != nL){ 
     #     p_val = wilcox.test(auc ~ risk.model, data = results_missForest_RF_AT, 
      #                                  paired = FALSE, conf.level = conf.level, alternative = "greater")[["p.value"]]
      #  }
      #}
      meanPercArea[unique(LABEL)==nL,nC] = mean(data[LABEL == nL,trueNC])
      stderrPercArea[unique(LABEL)==nL,nC] = sd(data[LABEL == nL,trueNC])/sqrt(sum(LABEL==nL))
      cat(colnames(data)[trueNC], ' (', nL , ') = ', meanPercArea[unique(LABEL)==nL,nC] , ' +- ', stderrPercArea[unique(LABEL)==nL,nC], '\n' ) 

    }
    
   
    
    
  
  }
  
  
}