p.value.analysis <- function(data, label,  var_types = NULL){
  source(file.path('.', 'Utils.R'))
  
  lab.levels = as.factor(label);
  if (is.null(var_types)){ 
    var_types = extract_type(data)
  }
  data.cat = data[c(grep('ord', var_types), grep('cat', var_types)),]
  
  data.num = data[-c(grep('ord', var_types), grep('cat', var_types)),]
  
  colnames.categorical = c("perc.on", paste("perc.on.in.", lab.levels[0], sep =""), paste("perc.on.in.", lab.levels[1], sep =""),
                           "p-value (chi square)", "direction")  
  df_categorical = as.data.frame(matrix(nrow = ncol(data.cat), ncol = length(colnames.categorical)))
    
  colnames.numeric = c("mean", "se", "min", "max", "p-value (<)", "p-value (>)", "direction")  
  df_numeric = as.data.frame(matrix(nrow = ncol(data.num), ncol = length(colnames.numeric)))
  for (nc in ncol(data)){
    var.now = data[,nc]
    idx_del = is.na(var.now)
    lab = label[-idx_del]
    var.now = var.now[-idx_del]
    
    if(strcmpi(unlist(var_types[nc]),'ord') || strcmpi(unlist(var_types[nc]),'cat')){
      p.value = chisq.test(x = var.now, y = lab, correct = TRUE)
      cat(colnames(data)[nc], " p-value = ", p.value[["p.value"]], '\n')
      perc_on.in.low = round(100*length(which(var.now ==1))/ length(var.now),1)
      perc_on.in.low = round(100*length(which(var.now ==1 & lab==lab.levels[0]))/ length(which(lab==lab.levels[0])),1)
      perc_on.in.high = round(100*length(which(var.now ==1 & lab==lab.levels[1]))/ length(which(lab==lab.levels[1])),1)
      cat(colnames(data)[nc], ';', perc_on.in.low, '(',length(which(var.now ==1 & lab==lab.levels[0])), ')', 
          ';', perc.on.in.high, '(',length(which(var.now ==1 & lab==lab.levels[1])), ')', '\n')
      
    }else{
      true = TRUE
      if (strcmpi(unlist(var_types[nc]),'num') ) {
        mean_v = round(mean(var.now, na.rm = true),2)
        se_v = round(sd(var.now, na.rm = true)/sqrt(length(which(is.finite(var.now)))),2)
        range_v = c(min(var.now, na.rm = true),max(var.now, na.rm = true)) 
        
        mean_low = round(mean(var.now[which(label ==lab.levels[0])], na.rm = true),2)
        se_low = round(sd(var.now[which(label ==lab.levels[0])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(label ==lab.levels[0])])))),2)
        range_low = c(min(var.now[which(label ==lab.levels[0])], na.rm = true),max(var.now[which(label ==lab.levels[0])], na.rm = true)) 
        
        mean_high = round(mean(var.now[which(label ==lab.levels[1])], na.rm = TRUE),2)
        se_high = round(sd(var.now[which(label ==lab.levels[1])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(label ==lab.levels[1])])))),2)
        range_high = c(min(var.now[which(label ==lab.levels[1])], na.rm = true),max(var.now[which(label ==lab.levels[1])], na.rm = true)) 
      }else{
        mean_v = median(var.now, na.rm = true)
        se_v = round(1.253*sd(var.now, na.rm = true)/sqrt(length(which(is.finite(var.now)))),2)
        range_v = c(min(var.now, na.rm = true),max(var.now, na.rm = true)) 
        
        mean_low = median(var.now[which(label ==lab.levels[0])], na.rm = true)
        se_low = round(1.253*sd(var.now[which(label ==lab.levels[0])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(label ==lab.levels[0])])))),2)
        range_low = c(min(var.now[which(label ==lab.levels[0])], na.rm = true),max(var.now[which(label ==lab.levels[0])], na.rm = true)) 
        
        mean_high = median(var.now[which(label ==lab.levels[1])], na.rm = true)
        se_high = round(1.253*sd(var.now[which(label ==lab.levels[1])], na.rm = true)/sqrt(length(which(is.finite(var.now[which(label ==lab.levels[1])])))),2)
        range_high = c(min(var.now[which(label ==lab.levels[1])], na.rm = true),max(var.now[which(label ==lab.levels[1])], na.rm = true)) 
        
      }
      cat(colnames(data)[nc], ';', mean_v, ' ±',  se_v, paste('[',range_v[1],', ',range_v[2],']', sep =""),';',
          mean_low, ' ±',  se_low, paste('[',range_low[1],', ',range_low[2],']', sep =""),';',
          mean_high, ' ±',  se_high, paste('[',range_high[1],', ',range_high[2],']', sep =""),'\n')
      
      pgreat= wilcox.test(var.now[label==lab.levels[0]],var.now[label==lab.levels[1]], alternative = 'greater', correct = TRUE)
      pless =  wilcox.test(var.now[label==lab.levels[0]],var.now[label==lab.levels[1]], alternative = 'less', correct = TRUE)
      if (pgreat$p.value < pless$p.value){ 
         cat(colnames(data)[nc], "class", lab.levels[0], "> class", lab.levels[1], pgreat$p.value, '\n')
      }else{
         cat(colnames(data)[nc], "class", lab.levels[0], "< class", lab.levels[1], pless$p.value, '\n')
      }
    }
    
  }

}

#  estraggo cutpoints per ogni punto

LABEL = as.factor(lab.use>=2)
data = data.use


results_with_radio = as.data.frame(matrix(0, ncol = 10, nrow =nrow(mean.imp)))
colnames( results_with_radio) = c('direction', 'cut-point','sensitivity', 'specificity', 'Acc','PPV', 'NPV','F1 score', 'AUC', 'MCC' )
results_with_radio$direction =  as.character(results_with_radio$direction)
n.boot.runs =10
for (nv in rownames(mean.imp)){
  cat(nv,':', '\n')
  var.now = data.use[,nv] #grep(nv, data.use)]
  
  if (is.ordered(var.now)){ chisq.test(x = var.now, y = LABEL, correct = TRUE)
  }else{
    pgreat= wilcox.test(var.now[lab.use==1],var.now[lab.use==0], alternative = 'greater', correct = TRUE)
    pless =  wilcox.test(var.now[lab.use==1],var.now[lab.use==0], alternative = 'less', correct = TRUE)
    if (pgreat$p.value < pless$p.value){ 
      direction = ">="
      results_with_radio[nv, 'direction'] = ">"
      cat("class 0 > class 1", pgreat$p.value, '\n')
    }else{
      direction = "<="
      results_with_radio[nv, 'direction'] = "<"
      
      cat("class 0 < class 1",pless$p.value, '\n')
    }
  }
  
  
  opt_cut <- cutpointr(x = var.now, class = as.factor(lab.use ==1),
                       silent = TRUE,
                       boot_runs = n.boot.runs,
                       boot_stratify = TRUE,
                       method = maximize_boot_metric,
                       metric = youden,
                       na.rm = TRUE, direction = direction)
  cutpoint = mean(opt_cut$boot[[1]]$optimal_cutpoint)
  
  
  TP = mean(opt_cut$boot[[1]]$TP_b)
  FP = mean(opt_cut$boot[[1]]$FP_b)
  TN = mean(opt_cut$boot[[1]]$TN_b)
  FN = mean(opt_cut$boot[[1]]$FN_b)
  sens =  mean(opt_cut$boot[[1]]$sensitivity_b)
  spec =  mean(opt_cut$boot[[1]]$specificity_b)
  Acc = (TP+TN)/(TP+FP+TN+FN)
  PPV = TP/(TP+FP)
  NPV = TN/(TN+FN)
  F1.score = 2*(PPV *sens)/(PPV+sens)
  MCC = ((TP*TN) - (FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  AUC =  mean(opt_cut$boot[[1]]$AUC_b)
  results_with_radio[nv, 2:ncol(results_with_radio)] =  
    c(cutpoint, sens, spec, Acc, PPV, NPV, F1.score, AUC, MCC)
  cat(direction,'\t', cutpoint, 
      '\t', sens,
      '\t', spec,
      '\t', Acc,
      '\t', PPV,
      '\t', NPV,
      '\t', F1.score,
      '\t', Acc,
      '\t',AUC,
      '\t',  MCC, '\n')
}



write.csv(learner, file=file.path(nameResults, paste('new_radio.score_maggiore7','.csv', sep="") ))

