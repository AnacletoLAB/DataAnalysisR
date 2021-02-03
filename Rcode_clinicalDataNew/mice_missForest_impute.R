mice_missForest_impute <- function(data, LABEL, method_impute = method_impute, dirData = file.path('.') ,
                                   num_imputations = 11, 
                                   visitOrder = "decreasing", maxiter = 20, ntree = 31){
  
  
  library(abind)
  library(mice)
  library(missForest)
  library(MissMech)
  library(caret)
  cat('Data Imputation', '\n')
  
  if (visitOrder == "decreasing"){
    decreasing = TRUE
    visitSequence = "revmonotone"
  }else{
    decreasing = FALSE
    visitSequence = "monotone"
  }
  
  LABEL = as.numeric(LABEL)-min(as.numeric(LABEL))
  
  if (length(grep("missForest", method_impute))>0){
  
    cat('\n', num_imputations, " Imputations with missForest",  '\n')
    imputed = vector("list", length = num_imputations)
    for (s in 1:num_imputations){
      nameDf = file.path(dirData, paste('data', strdata, '-imputed_missRF_', visitOrder, '_', as.character(s), '.csv', sep="",collapse = NULL))
      if (!file.exists(nameDf)){ 
        data_imputed_miss <- missForest(data, maxiter = maxiter, ntree = ntree,mtry = (ncol(data)-1)/3,
                                        variablewise = TRUE, decreasing =  decreasing)
        imp.num.n <- as.matrix(cbind("LABEL"=LABEL,data_imputed_miss$ximp))
        write.table(imp.num.n, file = nameDf, 
                    append = FALSE, quote = TRUE, sep = ",", 
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)  
      }else{
        imp.num.n <- as.matrix(read.table(file = nameDf, sep = ",", header = TRUE))
      }
      
      imputed[[s]] <-  imp.num.n
    }
    
  }
  
  
  if (length(grep("miceRF", method_impute))>0){
    # MICE imputation with random forest
    cat('\n', num_imputations, " Imputations with miceRF - maxiter = ", maxiter, '\n')
    nameDf = file.path(dirData, paste('data', strdata, '-imputed_miceRF_', visitOrder, '_', as.character(num_imputations),  '.csv', sep="",collapse = NULL))
    if (!file.exists(nameDf)){
      imp <- mice(data, m = num_imputations, method = "rf", maxit = maxiter,
                  ntree = ntree, mtry = (ncol(data)-1)/3,  visitSequence = visitSequence, printFLAG = FALSE)
    }else{ cat("miceRF: data already imputed!",  "\n", nameDf, "\n")}
    imputed = vector("list", length = num_imputations)
    for (s in 1:num_imputations){
      nameDf = file.path(dirData, paste('data', strdata, '-imputed_miceRF_', visitOrder, '_', as.character(s), '.csv', sep="",collapse = NULL))
      if (!file.exists(nameDf)){
        
        imp.num.n <- as.matrix(cbind("LABEL"=LABEL,complete(imp, action = s)))
        write.table(imp.num.n, file = nameDf, 
                    append = FALSE, quote = TRUE, sep = ",", 
                    eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
      }else{
         # cat("miceRF: read imputed data ",  "\n", nameDf, "\n")
          imp.num.n <- as.matrix(read.table(file = nameDf, sep = ",", header = TRUE))
      }
      imputed[[s]]  <-  imp.num.n
    }
  }
    
  if (length(grep("micePMM", method_impute))>0 | length(grep("miceNorm", method_impute))>0){
    cat('\n', num_imputations, " Imputations with micePMM - maxiter = ", maxiter, '\n')
    #imputo con parametri di default tranne maxiter : iù è altto più le impuatzioni convergono
    nameDf = file.path(dirData, paste('data', strdata, '-imputed_micePMM_', visitOrder, '_', as.character(num_imputations), '.csv', sep="",collapse = NULL))
    if (!file.exists(nameDf)){
      imp <- mice(data, m = num_imputations, method = "pmm", 
                  maxit = maxiter, visitSequence = visitSequence, printFLAG = FALSE)
    }else{ cat('\n', "already performed imputations:", "\n", nameDf, '\n')}
    
    imputed = vector("list", length = num_imputations)
    for (s in 1:num_imputations){
      nameDf = file.path(dirData, paste('data', strdata, '-imputed_micePMM_', visitOrder, '_', as.character(s), '.csv', sep="",collapse = NULL))
      if (!file.exists(nameDf)){
        #cat('\n', s, nameDf, '-', 'imputation not existent', '\n')
        imp.num.n <- ( cbind("LABEL" = LABEL,  as.matrix(complete(imp, action = s))))
        write.table(imp.num.n, file = nameDf, 
                    append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
      }else{
        #cat('\n', s, 'read file', nameDf, '\n')
        imp.num.n <- as.matrix(read.table(file = nameDf, sep = ",", header = TRUE))
      }
     
      imputed[[s]]  <-  imp.num.n
      
    }
  }   
    
  if (length(grep("distFree", method_impute))>0){
    cat('\n', num_imputations, " Imputations with distFree",  '\n')
    imputed = vector("list", length = num_imputations)
    #imputo con parametri di default tranne maxiter : iù è altto più le impuatzioni convergono
    for (s in 1:num_imputations){
      nameDf = file.path(dirData, paste('data', strdata, '-imputed_distFree_', visitOrder, '_', as.character(s), '.csv', sep="",collapse = NULL))
      if (!file.exists(nameDf)){
        no.na = t(apply(apply(apply(data,2, is.na),2,as.numeric),2,sum))
        imp <- Impute(apply(data[,order(no.na, decreasing = (visitOrder=="decreasing"))], 2, as.numeric), imputation.method = "dist.free")
        imp.num.n <- cbind("LABEL" = LABEL, as.matrix(imp$yimp))
        write.table(imp.num.n, file = nameDf,
                    append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
      }else{
        imp.num.n <- as.matrix(read.table(file = nameDf, sep = ",", header = TRUE))
      }
    
      imputed[[s]] <-  imp.num.n
      
      
    } 
    
  }
    
  return(imputed)
    
}