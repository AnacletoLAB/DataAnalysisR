p.value.analysis <- function( data, LABEL,  n.boot.runs = 100){    
  # prende il df dei dati contenenti un problema binario e da in outpput una serie di statistiche 
  # per le var categoriche e binarie clacola il p-value con chi square
  # per le variabili numeriche usa one-sided Mann Whitney U test (il p-value è il minore tra il one sided < e one-sided per il >)
  # LABEL è la colonna delle LABEL
  # n.boot.runs = numero di bootstrapped samples per stimare l'optimal cutpoint; se  n.boot.runs =1 
  # l'optimal cutpoint è stimato su tutto il campione
  # il tipo delle features viene inferito dal nome:
  # nomi che cominciano con :
  # INT.* : sono variabili a valori interi
  # NUM.* : variabili a valori reali
  # BIN.* : sono variabili a valori categorici BINARI 
  # CAT.* : sono variabili a valori categorici  
  # ORD.* : sono variabili a valori ordinali  
  
 print("Statistical analysis of input data \n")
  var_types <- extract_type(data)
  lab.levels = levels(LABEL)
  if (length(lab.levels)>2){
    print('only binary levels are analyzed')
  }
  # 1 row per ogni feature
  
  
  
 
  
  
  colnames.df = c("all", paste("class ", lab.levels[1], sep =""), paste("class ", lab.levels[2], sep =""), 
                       "p-value",
                       "direction", "sens", "spec", "acc","PPV", "NPV", "MCC", "F1", "Auc")
  df = as.data.frame(matrix(nrow = ncol(data), ncol = length(colnames.df)))
  rownames(df ) = colnames(data)
  colnames(df) = colnames.df
  
    
  for (nc in 1:ncol(data)){
      var.now = data[,nc]
      idx_use = !is.na(var.now)
      lab = LABEL[idx_use]
      var.now = var.now[idx_use]
      var.now_numeric= var.now
      if(pracma::strcmp(unlist(var_types[nc]),'bin') || pracma::strcmp(unlist(var_types[nc]),'ord') || pracma::strcmp(unlist(var_types[nc]),'cat')){
        var.now = as.factor(var.now)
        levels_var = levels(as.factor(var.now))
        label_on = levels_var[length(levels_var)]
        cat('for categorical variable on value is = ', label_on, '\n')
        
        p.value = chisq.test(x = var.now, y = lab, simulate.p.value = TRUE)
        
        str_chi = paste("p-value = ", round(p.value[["p.value"]],4), '\n', sep ="")
        
        perc.on = round(100*length(which(var.now == label_on))/ length(var.now),1)
        perc.on.in.low = round(100*(length(which(var.now == label_on & lab==lab.levels[1]))/ length(which(lab==lab.levels[1]))),1)
        perc.on.in.high = round(100*(length(which(var.now == label_on & lab==lab.levels[2]))/ length(which(lab==lab.levels[2]))),1)
        # cat(colnames(data)[nc], '\t'  , 
        #     perc.on,'(', length(which(var.now == label_on)),')', '\t'  ,
        #     perc.on.in.low, '(',length(which(var.now ==label_on & lab==lab.levels[1])), ')', '\t'  ,
        #     perc.on.in.high, '(',length(which(var.now ==label_on & lab==lab.levels[2])), ')', '\t', str_chi)
        if (perc.on.in.high>=perc.on.in.low){ 
          direction =paste("class ", lab.levels[2], " >= class ", lab.levels[1])
          est_dir = ">="
        }else{
          est_dir = "<="
          direction = paste("class ", lab.levels[1], " >= class ", lab.levels[2])
        }
        first_estimates  = c(paste(perc.on, '%', sep=""), paste(perc.on.in.low, '%', sep=""),  paste(perc.on.in.high, '%', sep=""), round(p.value[["p.value"]],4), direction)
        
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
        
        
        pgreat= wilcox.test(data[lab==lab.levels[1],nc],data[lab==lab.levels[2],nc],na.action = "na.omit", alternative = 'greater',  PAIRED = FALSE)
        pless =  wilcox.test(data[lab==lab.levels[1],nc],data[lab==lab.levels[2],nc], na.action = "na.omit", alternative = 'less',  PAIRED = FALSE)
        if (pgreat$p.value < pless$p.value){ 
          str_w = paste("class ", lab.levels[1], " > class ", lab.levels[2], " p-value = ", round(pgreat$p.value,4), '\n', sep ="")
          direction = paste("class ", lab.levels[1], " > class ", lab.levels[2])
          est_dir ="<="    # per cutpointr x is supposed to be smaller for the positive class
          p.val = pgreat$p.value
        }else{
          est_dir = ">="   # per cutpointr x is supposed to be smaller for the positive class
          str_w = paste("class ", lab.levels[1], " < class ", lab.levels[2], " p-value = ", round(pless$p.value,4), '\n', sep ="")
          direction = paste("class ", lab.levels[1], " > class ", lab.levels[2])
          p.val = pless$p.value
        }
        
       
      
        first_estimates = c(paste( mean_v, ' ± ',  se_v, ' [',range_v[1],', ',range_v[2],']', sep =""),  
          paste(mean_low, ' ± ',  se_low, ' [',range_low[1],', ',range_low[2],']', sep =""),
          paste(mean_high, ' ± ',  se_high,' [',range_high[1],', ',range_high[2],']', sep =""), 
          round(p.val,4),  direction)
      }
      
      if (n.boot.runs==1){ 
        boot_stratify = FALSE
      }else{ 
        boot_stratify = TRUE
      }
      
      opt_cut <- cutpointr::cutpointr(x = var.now_numeric, class = as.factor(lab ==1),
                           silent = TRUE,
                           boot_runs = n.boot.runs,
                           boot_stratify = boot_stratify,
                           method = maximize_boot_metric,
                           metric = youden,
                           na.rm = TRUE, direction = est_dir)
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
     
      # cat(direction,'\t', cutpoint, 
      #     '\t', sens,
      #     '\t', spec,
      #     '\t', Acc,
      #     '\t', PPV,
      #     '\t', NPV,
      #     '\t', F1.score,
      #     '\t', Acc,
      #     '\t',AUC,
      #     '\t',  MCC, '\n')
      # 
      first_estimates = c(first_estimates, sens,spec, Acc,PPV,NPV, MCC, F1.score, AUC)
      print(first_estimates)
      
      df[colnames(data)[nc], ] = first_estimates
      
     
      
      
  }
  
  return(df)
  
}