check_MCAR <- function(data, LABEL,  nIt = 50 , thrVar = 0.025) {
  
  library(MissMech)
  library(missForest)
  library(VIM)
  library(pracma)
  source(file.path('.','RF_on_RF_and_Boruta.R'))
  #thrVar = 0.025 # variable whose covariance is lower than this value are removed for they have to low variance
  #nIt = 50 #numero di imputazioni per missForest
  
  
  cat('Data Imputation', '\n')
  
  #srimuovo i sample con troppi NAN
  sumNA <-  apply(is.na(data), 1, sum)
  remove_samples_too_nan <-  which(sumNA/ncol(data)>1/2) 

  if (length(remove_samples_too_nan)>0){
    cat('Removed for too many missing values: ', '\n',rownames(data)[remove_samples_too_nan], '\n')
    data <- data[-remove_samples_too_nan]
  }
  
    
  #rimuovo le features con troppi NAN
  sumNA <-  apply(is.na(data), 2, sum)
  remove_too_nan <-  which(sumNA/nrow(data)>1/2) 

  if (length(remove_too_nan)>0){
    cat('Removed for too many missing values: ', '\n',colnames(data)[remove_too_nan], '\n')
    data <- data[,-remove_too_nan]
  }
  
  # calcolo la covarianza di ogni variabile e la rimuovo se ha troppa poca variabilita
  resCov = round(cov(apply(data,2,as.numeric), use = "pairwise.complete.obs"),3)
  # la diagonale mi da la covarianza della variabile
  idxDel = which(diag(resCov)<thrVar)
  if (length(idxDel)>0){
    cat('Removed for too low variability: ',  '\n', colnames(data)[idxDel], '\n')
    data = data[, -idxDel]
  }
  
  
  for (nc in 1: length(levels(LABEL))){
    cat("class ", levels(LABEL)[nc], " has ", length(which(LABEL == levels(LABEL)[nc])),"samples\n")
  }

  
  
  if (any(is.na(data))){
    # imputo con missForest in increasing order: mi servono per passarli a check_MCAR perchè
    # speddo distFree si impalla
   
    #innanzi tutto valuto se sia meglio increasing o decreasing order, in base a OOB error
    OOB.inc  = 0
    OOB.dec  = 0
    imp.inc = vector("list", length = nIt)
    sum.inc = matrix(0,nrow = nrow(data), ncol = ncol(data))
    imp.dec = vector("list", length = nIt)
    sum.dec = matrix(0,nrow = nrow(data), ncol = ncol(data))
    
    for (it in 1:nIt){
      data_imputed_miss.inc <- missForest(xmis=data, variablewise = TRUE, decreasing =  FALSE)
      OOB.inc = OOB.inc+mean(data_imputed_miss.inc$OOBerror)  
      imp.inc[[it]] <- cbind("LABEL"=LABEL,data_imputed_miss.inc$ximp)
      sum.inc = sum.inc+data_imputed_miss.inc$ximp
      data_imputed_miss.dec <- missForest(xmis=data, variablewise = TRUE, decreasing =  TRUE)
      OOB.dec = OOB.dec+mean(data_imputed_miss.dec$OOBerror)
      imp.dec[[it]] <- cbind("LABEL"=LABEL,data_imputed_miss.dec$ximp)
      sum.dec = sum.dec+data_imputed_miss.dec$ximp
    }
    data.imp = data
    
    if(OOB.inc <OOB.dec) { 
      print('Increasing imputation order used\n')
      data.imp[is.na(data)] = sum.inc[is.na(data)]/nIt 
    }else{
      print('Decreasing imputation order used\n')
      data.imp[is.na(data)] = sum.dec[is.na(data)]/nIt 
    }
    
    data.imp = apply(transform.data.with.type(data.imp),2,as.numeric)
    
    
    res = tryCatch(TestMCARNormality(data = dataMCAR, del.lesscases = 1, imputed.data =data.imp ), 
                   error = function(e){print(e)},
                   finally = print("data is MCAR"))
  }else{
      data.imp = data
      }  
    
  dataMCAR = data
  names(dataMCAR) = colnames(data)
  par(pin=c(0,0))
  data_aggr = VIM::aggr(dataMCAR, col=c("skyblue", "red", "orange"), numbers=TRUE, sortVars=TRUE, 
                        cex.axis=.7, 
                        ylab=c("Proportion of missingness","Missingness Pattern"))
    
  
  
  return(data.imp)

}