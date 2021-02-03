check_MCAR <- function(data, LABEL,  nIt = 50 , thrVar = 0.025, 
                       imputed_data= imp) {
  
  library(MissMech)
  library(VIM)
  library(pracma)
  source(file.path('.','RF_on_RF_and_Boruta.R'))
  
  
  
  cat('Data Imputation', '\n')
  
  
  
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
    res = tryCatch(TestMCARNormality(data = data, del.lesscases = 1, imputed.data =imputed_data ), 
                   error = function(e){print(e)},
                   finally = print("data is MCAR"))
  }else{
      imputed_data = data
  }  
    
  dataMCAR = data
  names(dataMCAR) = colnames(data)
  par(pin=c(0,0))
  data_aggr = VIM::aggr(dataMCAR, col=c("skyblue", "red", "orange"), numbers=TRUE, sortVars=TRUE, 
                        cex.axis=.7, 
                        ylab=c("Proportion of missingness","Missingness Pattern"))
    
  
  
  return(res)

}