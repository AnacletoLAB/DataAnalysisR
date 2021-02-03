plot.variances.imputation.methods <- function(data = data, LABEL = LABEL, ntree = 31, maxiter =101, num_imputs =101, strdata='data', dirData = 'data') {
  
  
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
  
  method_imputation = c( "missForest", "miceRF", "micePMM")
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
      
      result.df[(grep(method_impute, method_imputation)-1)*length(visit_orders)+n_order, ] = 
        c(method_impute, visit_orders[n_order],  round(mean(variances_el),prec_digit),
          round(sd(variances_el)/sqrt(length(variances_el)),prec_digit), 
          round(min(variances_el),prec_digit), 
            round(max(variances_el),prec_digit))
    }  
        
    if (grep(method_impute, method_imputation)==1){
      mega_df = df
    }else{
      mega_df = rbind(mega_df, df)
    }
      

  }
  
  colnames(mega_df) = c("imputation.method", "no.imputations", "increasing", "decreasing")
  ggplot() +  
    geom_line(data= mega_df, aes(x=no.imputations, y=increasing), linetype = "dashed")+
    geom_point(data= mega_df, aes(x = no.imputations, y=increasing, color = increasing)) +
    geom_line(data= mega_df, aes(x=no.imputations, y=decreasing), linetype = 4)+
    geom_point(data= mega_df, aes(x = no.imputations, y=decreasing, color = decreasing), shape = 25)+ 
    ylab("global Variance")+xlab("no. imputations")+
    facet_wrap(~imputation.method)+theme_bw()+ 
    labs(
      shape = "boh"
    )
  
  p = vector("list", length = length(method_imputation))
  
  p[[1]] <- ggplot() +  
    geom_line(data= mega_df[which(mega_df$imputation.method=="miceRF"),], aes(x=no.imputations, y=increasing), linetype = "dashed")+
    geom_point(data= mega_df[which(mega_df$imputation.method=="miceRF"),], aes(x = no.imputations, y=increasing, color = increasing))+
    geom_line(data= mega_df[which(mega_df$imputation.method=="miceRF"),], aes(x=no.imputations, y=decreasing), linetype = 4)+
    geom_point(data= mega_df[which(mega_df$imputation.method=="miceRF"),], aes(x = no.imputations, y=decreasing, color = decreasing), shape = 25)+
    ylab("global var")+ggtitle("miceRF")
  
  p[[2]] <- ggplot() +  
    geom_line(data= mega_df[which(mega_df$imputation.method=="micePMM"),], aes(x=no.imputations, y=increasing), linetype = "dashed")+
    geom_point(data= mega_df[which(mega_df$imputation.method=="micePMM"),], aes(x = no.imputations, y=increasing, color = increasing))+
    geom_line(data= mega_df[which(mega_df$imputation.method=="micePMM"),], aes(x=no.imputations, y=decreasing), linetype = 4)+
    geom_point(data= mega_df[which(mega_df$imputation.method=="micePMM"),], aes(x = no.imputations, y=decreasing, color = decreasing), shape = 25)+
    ylab("global var")+ggtitle("micePMM")
  
  p[[3]] <- ggplot() +  
    geom_line(data= mega_df[which(mega_df$imputation.method=="missForest"),], aes(x=no.imputations, y=increasing), linetype = "dashed")+
    geom_point(data= mega_df[which(mega_df$imputation.method=="missForest"),], aes(x = no.imputations, y=increasing, color = increasing))+  
    geom_line(data= mega_df[which(mega_df$imputation.method=="missForest"),], aes(x=no.imputations, y=decreasing), linetype = 4)+
    geom_point(data= mega_df[which(mega_df$imputation.method=="missForest"),], aes(x = no.imputations, y=decreasing, color = decreasing), shape = 25)+
    ylab("global var")+ggtitle("missForest")
  
  grid.arrange(p[[1]], p[[2]], p[[3]], nrow = 2
               )
  
  
  ggplot() +  
    geom_line(data= mega_df, aes(x=no.imputations, y=increasing), linetype = "dashed")+
    geom_point(data= mega_df, aes(x = no.imputations, y=increasing, color = increasing)) +
    geom_line(data= mega_df, aes(x=no.imputations, y=decreasing), linetype = "dotdash")+
    geom_point(data= mega_df, aes(x = no.imputations, y=decreasing, color = decreasing))+ 
    ylab("global Variance")+xlab("no. imputations")+
    facet_wrap(~imputation.method)
  
  
  
  # ora plotto!
  myPalette = brewer.pal(n = 9, name = "Blues")
  p = vector("list", length = length(method_imputation))
    
  miny_mice = min(c(variances[[grep("miceRF",method_imputation)]], variances[[grep("micePMM",method_imputation)]])) 
  maxy_mice = max(c(variances[[grep("miceRF",method_imputation)]], variances[[grep("micePMM",method_imputation)]])) 
  rng = c(min(unlist(variances)), max(unlist(variances)))
  
  for (no_m_impute in 1:length(method_imputation)){
    method_impute = method_imputation[no_m_impute]
    mean_vars = variances[[no_m_impute]]
    df = data.frame("noImputations" = 2:num_imputs, "meanVar" = mean_vars)
    if (method_impute == "micePMM"){
      tit = paste("micePMM", " (", visit_orders[j], " order)", sep ="")
    }else{
      tit = paste(method_impute, " (", visit_orders[j], " order)", sep ="")
    }
    
    if (method_impute=="miceRF" || method_impute=="micePMM"){     

      p[[grep(method_impute, method_imputation)]] <- ggplot(data= df, aes(x=noImputations, y=meanVar, group=1)) +
      geom_line(linetype = "dashed")+ylim(miny_mice-(maxy_mice-miny_mice)/10, maxy_mice+(maxy_mice-miny_mice)/10)+
      geom_point(aes(color = meanVar))+ ggtitle(tit)+ylab("mean Variance")+xlab("no. imputations")
    }else{
      p[[grep(method_impute, method_imputation)]] <- ggplot(data=df, aes(x=noImputations, y=meanVar, group=1)) +
        geom_line(linetype = "dashed")+
        geom_point(aes(color = meanVar))+ ggtitle(tit)+ylab("mean Variance")+xlab("no. imputations")
      
    }
  }
  grid.arrange(arrangeGrob(  p[[1]], p[[2]],p[[3]],  nrow = 2, ncol = 2))

  return(list(mega_df, result_df))
 
 
}  #endfun




