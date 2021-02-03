check_MCAR <- function(data, LABEL, thrVar = 0.025) {
  
  library(MissMech)
  library(VIM)
  library(pracma)
  source(file.path('.','RF_on_RF_and_Boruta.R'))
 
  cat('Data Imputation', '\n')
  sumNA <-  apply(is.na(data), 1, sum)
  remove_samples_too_nan <-  which(sumNA/ncol(data)>1/2) 
  
  sumNA <-  apply(is.na(data), 2, sum)
  remove_too_nan <-  which(sumNA/nrow(data)>1/2) 
  # 
  
  if (length(remove_too_nan)>0){
    cat('Removed for too many missing values: ', '\n',colnames(data)[remove_too_nan], '\n')
    data <- data[,-remove_too_nan]
  }
  
  resCov = round(cov(apply(data,2,as.numeric), use = "pairwise.complete.obs"),3)
  idxDel = which(diag(resCov)<thrVar)
  if (length(idxDel)>0){
    cat('Removed for too low variability: ',  '\n', colnames(data)[idxDel], '\n')
    data = data[, -idxDel]
  }
  
  sumNA_row <-  round(apply(is.na(data), 1, sum)/ncol(data),2)
  del_point = which(sumNA_row>0.25)
  if (length(del_point)>0){
    cat(length(del_point) , " samples removed for they had more than 25% of missing values")
    data = data[-del_point,]
    LABEL = LABEL[-del_point] 
  }
  
  for (nc in 1: length(levels(LABEL))){
    cat("class ", levels(LABEL)[nc], " has ", length(which(LABEL == levels(LABEL)[nc])),"samples")
  }

  dataMCAR = apply(data, 2, as.numeric)
  
  res = TestMCARNormality(data = dataMCAR, imputation.method = "Dist.Free")
  
  names(dataMCAR) = colnames(data)
  par(pin=c(0,0))
  data_aggr = VIM::aggr(dataMCAR, col=c("skyblue", "red", "orange"), numbers=TRUE, sortVars=TRUE, 
                      cex.axis=.7, 
                      ylab=c("Proportion of missingness","Missingness Pattern"))
  

  var_types <- extract_type(data)
  lab.levels = levels(LABEL)
  for (nc in 1:ncol(data)){
    var.now = data[,nc]
    idx_use = !is.na(var.now)
    lab = LABEL[idx_use]
    var.now = var.now[idx_use]
    
    if(strcmp(unlist(var_types[nc]),'ord') || strcmp(unlist(var_types[nc]),'cat')){
      p.value = chisq.test(x = var.now, y = lab, simulate.p.value = TRUE)
      str_chi = paste("p-value = ", round(p.value[["p.value"]],4), '\n', sep ="")
      perc.on = round(100*length(which(var.now ==1))/ length(var.now),1)
      perc.on.in.low = round(100*(length(which(var.now ==1 & lab==levels(lab)[1]))/ length(which(lab==lab.levels[1]))),1)
      perc.on.in.high = round(100*(length(which(var.now ==1 & lab==levels(lab)[2]))/ length(which(lab==lab.levels[2]))),1)
      cat(colnames(data)[nc], '\t'  , 
          perc.on,'(', length(which(var.now ==1)),')', '\t'  ,
          perc.on.in.low, '(',length(which(var.now ==1 & lab==lab.levels[1])), ')', '\t'  ,
          perc.on.in.high, '(',length(which(var.now ==1 & lab==lab.levels[2])), ')', '\t', str_chi)
      
    }else{
      true = TRUE
      if (strcmp(unlist(var_types[nc]),'num') ) {
        mean_v = round(mean(var.now, na.rm = true),2)
        se_v = round(sd(var.now, na.rm = true)/sqrt(length(which(is.finite(var.now)))),2)
        range_v = c(min(var.now, na.rm = true),max(var.now, na.rm = true)) 
        
        mean_low = round(mean(var.now[which(lab ==lab.levels[1])], na.rm = true),2)
        se_low = round(sd(var.now[which(lab ==lab.levels[1])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(lab ==lab.levels[1])])))),2)
        range_low = c(min(var.now[which(lab ==lab.levels[1])], na.rm = true),max(var.now[which(lab ==lab.levels[1])], na.rm = true)) 
        
        mean_high = round(mean(var.now[which(lab ==lab.levels[2])], na.rm = TRUE),2)
        se_high = round(sd(var.now[which(lab ==lab.levels[2])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(lab ==lab.levels[2])])))),2)
        range_high = c(min(var.now[which(lab ==lab.levels[2])], na.rm = true),max(var.now[which(lab ==lab.levels[2])], na.rm = true)) 
      }else{
        mean_v = median(var.now, na.rm = true)
        se_v = round(1.253*sd(var.now, na.rm = true)/sqrt(length(which(is.finite(var.now)))),2)
        range_v = c(min(var.now, na.rm = true),max(var.now, na.rm = true)) 
        
        mean_low = median(var.now[which(lab ==lab.levels[1])], na.rm = true)
        se_low = round(1.253*sd(var.now[which(lab ==lab.levels[1])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(lab ==lab.levels[1])])))),2)
        range_low = c(min(var.now[which(lab ==lab.levels[1])], na.rm = true),max(var.now[which(lab ==lab.levels[1])], na.rm = true)) 
        
        mean_high = median(var.now[which(lab ==lab.levels[2])], na.rm = true)
        se_high = round(1.253*sd(var.now[which(lab ==lab.levels[2])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(lab ==lab.levels[2])])))),2)
        range_high = c(min(var.now[which(lab ==lab.levels[2])], na.rm = true),max(var.now[which(lab ==lab.levels[2])], na.rm = true)) 
        
      }
     
      
      pgreat= wilcox.test(data[LABEL==lab.levels[1],nc],data[LABEL==lab.levels[2],nc],na.action = "na.omit", alternative = 'greater', correct = TRUE)
      pless =  wilcox.test(data[LABEL==lab.levels[1],nc],data[LABEL==lab.levels[2],nc], na.action = "na.omit", alternative = 'less', correct = TRUE)
      if (pgreat$p.value < pless$p.value){ 
        str_w = paste("class ", lab.levels[1], " > class ", lab.levels[2], " p-value = ", round(pgreat$p.value,4), '\n', sep ="")
      }else{
        str_w = paste("class ", lab.levels[1], " < class ", lab.levels[2], " p-value = ", round(pless$p.value,4), '\n', sep ="")
      }
      
      cat(colnames(data)[nc], '\t', mean_v, ' ±',  se_v, paste('[',range_v[1],', ',range_v[2],']', sep =""),'\t',
          mean_low, ' ±',  se_low, paste('[',range_low[1],', ',range_low[2],']', sep =""),'\t',
          mean_high, ' ±',  se_high, paste('[',range_high[1],', ',range_high[2],']', sep =""), '\t', str_w)
    }
    
  }
  

  
  
  return(data)

}