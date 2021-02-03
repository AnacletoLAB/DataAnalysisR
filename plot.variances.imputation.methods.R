plot.variances.imputation.methods <- function(data = data, LABEL = LABEL, ntree = 31, maxiter =101, num_imputs =101, strdata='data', dirData = 'data') {
  
 # restituisce un array con due elementi dal nome method_order: in method_order[1] ho il metodo di imputazione migliore
  # method_order[2] ho l'ordine di imputazione migliore
  
   conf.level = 0.05
  # imputo num_imputs volte i dati e calcolo la varianza sui dati imputati; scelgo il metodo che ha varianza minima
  library(pROC)
  library(caret)
  library(gridExtra)
  library(gtable)
  library(RColorBrewer)
  
  
 
  source(file.path('.', 'mice_missForest_impute.R'))
  
  source(file.path('.', 'Utils_preprocess.R'))
  source(file.path('.', 'Utils_postprocess.R'))
  source(file.path('.', 'rubins_rule.R'))

 
  nameLab = "LABEL"
  
  dirImpute = file.path(dirData,'Imputed') 
  if (!dir.exists(dirImpute)){
    dir.create(dirImpute)
  }
  
  method_imputation = c( "missRF", "miceRF", "micePMM")
  method_imputation_categorical = as.factor(method_imputation)
  
  # uso dati normalizzati perche' voglio dare una misura della varianza delle imputazioni!
  list_norm = min_max_norm(data)
  data_norm = list_norm[["norm"]]
  # mi tengo i coefficienti della normalizzazione per fare lo stesso sui dati imputati!
  col_ceff = list("min"=list_norm[["min"]],"max"=list_norm[["max"]]) 
  idx_nan = which(is.na(data))
 
 # computo le variance delle imputazioni usando da 5 a 50 imputazioni!
  visit_orders = c("increasing", "decreasing")
  mega_df = NULL
  prec_digit = 5
  imputed = vector("list", length=length(visit_orders))
  
  result.df = data.frame(matrix(ncol = 6, nrow = length(method_imputation)*length(visit_orders)))
  colnames(result.df) = c("method.imputation", "order", "variance", "se", "min", "max")
  variances = vector("list", length = length(method_imputation)*length(visit_orders))           
  for (method_impute in method_imputation_categorical){
    df = data.frame(imputation.method = rep(method_impute, length=num_imputs-2+1), 
                    no.imputations = 2:num_imputs)
        
    for (n_order in 1:length(visit_orders)){
      mini_df = data.frame(imputation.method = rep(method_impute, length=num_imputs-2+1), 
                           imputation.order = rep(visit_orders[n_order], length=num_imputs-2+1),
                           no.imputations = 2:num_imputs)
     # dirImpute = file.path(dirData, paste('Imputed', strdata,'_', visit_orders[n_order],'_', ntree, sep =""))
      imputed = mice_missForest_impute(method_impute = method_impute, data = data, LABEL = LABEL, 
                                     dirData = dirImpute, 
                                     maxiter = maxiter, num_imputations = num_imputs, ntree = ntree,  
                                     visitOrder = visit_orders[n_order], strdata = strdata)
    
      imputed_values = matrix(NA, nrow = length(idx_nan), ncol = num_imputs)
      variances_el = 0
      for (n_imp in 1:num_imputs){
        imputed_normalized = norm_with_coeffs(imputed[[n_imp]][,2:ncol(imputed[[1]])], col_ceff)
        imputed_values[,n_imp] = imputed_normalized[idx_nan]
        
        if (n_imp>1){
          variances_el = c(variances_el, mean(apply(imputed_values[,1:n_imp],1,var)))
        }
      }  
      
      variances_el = variances_el[2:length(variances_el)]
      variances[[(grep(method_impute, method_imputation)-1)*length(visit_orders)+n_order]] = variances_el
      cat(paste(method_impute, "\t", visit_orders[n_order], "\t",
          round(mean(variances_el),prec_digit), " +- ", round(sd(variances_el)/sqrt(length(variances_el)),prec_digit),
          " [", round(min(variances_el),prec_digit), " , ",
          round(max(variances_el),prec_digit),"]",
          "\n", sep=""))
      c.names = colnames(df)
      df = cbind(df,  variances_el)
      colnames(df) = c(c.names, visit_orders[n_order])
      #salvo il risultato nella colonna dei risultati
      c.names= colnames(mini_df)
      mini_df = cbind(mini_df, variances_el)
      colnames(mini_df) = c(c.names, "variance")
      
      result.df[(grep(method_impute, method_imputation)-1)*length(visit_orders)+n_order, ] = 
        c(method_impute, visit_orders[n_order],  round(mean(variances_el),prec_digit),
          round(sd(variances_el)/sqrt(length(variances_el)),prec_digit), 
          round(min(variances_el),prec_digit), 
            round(max(variances_el),prec_digit))
      
      if ((grep(method_impute, method_imputation)==1) & (n_order==1)){
        mega_df = mini_df
      }else{
        mega_df = rbind(mega_df, mini_df)
      }
    }  
        
   
      

  }
  
  pd <- position_dodge(width = 0.4)
  colnames(mega_df) = c("imputation.method",  "order", "no.imputations", "variance")
  gg <- ggplot(mega_df, aes=(group= imputation.method)) +  
    geom_line(data= mega_df, aes(x=no.imputations, y=variance, group= imputation.method))+
    geom_point(data= mega_df, aes(x = no.imputations, y=variance, color = order), size=0.1) +
#    geom_point(data= mega_df, aes(x = no.imputations, y=decreasing, color = decreasing), shape = 25)+ 
    ylab("global Variance")+xlab("no. imputations")+
    facet_wrap(~imputation.method)+theme_bw()

  ggsave(filename =file.path(dirData, 'compareVariances.jpg') , plot =gg)
  p = vector("list", length = length(method_imputation))
 
  for (method_impute in method_imputation_categorical){
  
    p[[grep(method_impute, method_imputation_categorical)]] <- ggplot(mega_df[which(mega_df$imputation.method==method_impute),], aes(x=no.imputations, y=variance, group = order)) + 
      geom_line(aes(color = order), size = .3, position = pd) +
      ggtitle(method_impute)+theme_bw()
     # theme(legend.position="bottom")
  }
  
 
  ml <-marrangeGrob(p, nrow=2, ncol=2)
  ml
  ggsave(filename =file.path(dirData, 'compareVariances_differentscale.jpg'), plot=ml)
  
  
  # wilcoxon test di ogni metodo + ordine di imputazione contro ogni altro metodo + ordine di imputazione
  # per vedere quale è minore degli altri
  wilcoxon_less = matrix(1, nrow = length(method_imputation)*length(visit_orders), ncol = length(method_imputation)*length(visit_orders))
  
  for (n_method_0 in 1:length(method_imputation_categorical)){
    method_impute_0 = method_imputation_categorical[n_method_0]  # metodo che confronteò con tutti gli altri
    for (n_order_0 in 1:length(visit_orders)){
      order_impute_0 = visit_orders[n_order_0]
      if (n_method_0==1 & n_order_0==1){
        col_row_names = paste(method_impute_0,  order_impute_0, sep="_")
      }else{
        col_row_names = c(col_row_names, paste(method_impute_0,  order_impute_0, sep="_"))
      }
      vars_0 = mega_df[which((mega_df$imputation.method==method_impute_0)  & (mega_df$order==order_impute_0)), ]   
      
      
      for (n_method_1 in 1:length(method_imputation_categorical)){
        method_impute_1 = method_imputation_categorical[n_method_1]
        for (n_order_1 in 1:length(visit_orders)){
          if ((n_method_1 != n_method_0) || (n_order_0 != n_order_1)){
            order_impute_1 = visit_orders[n_order_1]
            
            vars_1 = mega_df[which((mega_df$imputation.method==method_impute_1)  & (mega_df$order==order_impute_1)), ]   
            pvalL = wilcox.test(vars_0$variance, vars_1$variance, paired = TRUE, alternative = "less")
           
            wilcoxon_less[(n_method_0-1)*length(visit_orders)+n_order_0,(n_method_1-1)*length(visit_orders)+n_order_1 ] =pvalL$p.value

          }          
            
        }
      }
    }
  }
    
  colnames(wilcoxon_less)= col_row_names
  rownames(wilcoxon_less)= col_row_names
  
  # trovo il metodo che ha between imputation variance minore degli altri
  sum_less = apply(wilcoxon_less<conf.level, 1, sum)
  #best method and order are those in the neame of the column
  method_order = col_row_names[which(sum_less == max(sum_less))]
  method_order = strsplit(method_order, "_")[[1]]
  
  #in method_order[[1]][1] ho il metodo di imputazione migliore
  # method_order[[1]][2] ho l'ordine di imputazione migliore
 
  return(method_order)
 
}  #endfun




