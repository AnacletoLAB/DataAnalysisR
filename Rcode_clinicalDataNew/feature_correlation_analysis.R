feature_correlation_analysis <- function(nIter = 5, ntree = 500, num_folds = 10, num.inner.folds = 7, num.boruta.folds = 5, strDir = '_covnet_score',  
                                    strdata = '_covnet_score') {
  
  library(pROC)
  library(caret)
  library(gridExtra)
  library(gtable)
  library(RColorBrewer)
  library(ggpubr)
  
  dirSw = paste('C:/DATI/Policlinico-Humanitas-COVIRT/Andrea_RXCovid/MATLAB_CODE/data', strDir, '/Rcode_clinicalDataNew/', sep ="")
  setwd(dirSw)
  
  source(file.path('.', 'mice_missForest_impute.R'))
  source(file.path('.', 'RF_on_RF_and_Boruta.R'))
  source(file.path('.', 'glmnet_regression.R'))
  source(file.path('.','perf.meas.R'))
  
  source(file.path('.', 'Utils_preprocess.R'))
  source(file.path('.', 'Utils_postprocess.R'))
  source(file.path('.', 'rubins_rule.R'))
  
  source(file.path('.', 'check_MCAR.R'))
  source(file.path('.', 'find.opt.forest.R'))
  
  dirData = file.path('..')
  
  
  nIter = 3 # for each imputation run the process with 5 different 10-cv folds
  ntree = 101 # maximum number of trees used for data imputation
  num_folds = 10
  num.inner.folds = 5
  num.boruta.folds = 3
  maxiter = 11 # maximum number of imputation iterations
  visitOrder = "increasing"
  strDir = '_covnet_score'
  strdata = '_covnet_score'
  num_imputations = 25
  thrConf = 0.05
  
  #cambiare con nameLab = name_labels[1] per provare con  classi binarie
  nameLab = "LABEL"
  
  nameResults = file.path(dirData, paste('Results', strdata, '_', visitOrder, '_m_', num_imputations, sep =""))
  dirImpute = file.path(dirData, paste('Imputed', strdata,'_', visitOrder, '_', ntree, sep =""))
  data = read.csv(file = file.path(dirData, 
                                     paste('data', strdata, '.csv', sep="",collapse = NULL)),  
                    header = TRUE, sep = ",")
  
  
  if (!file.exists(nameResults)){
    dir.create(nameResults)
  }
  
  
  if (!file.exists(dirImpute)){
    dir.create(dirImpute)
  }
  
  LABEL = as.factor(data$LABEL)
  data = data[,-grep("LABEL", colnames(data))]
  metrics <- c("Sensitivity", "Specificity", "Pos Pred Value", "Neg Pred Value", "F1", "Balanced Accuracy" ) 
  folds.list = vector("list",length = nIter)
  for ( n.it in 1:nIter){
    ###########
    folds.list[[n.it]] <- caret::createFolds(y= LABEL, k = num_folds, list = TRUE, returnTrain = FALSE)
    ##########
  }
  
  data = check_MCAR(data = data, LABEL = LABEL, thrVar = 0.025)
  
  method_imputation = c(  "missForest")

  
 
  for (method_impute in method_imputation){
    cat("start with impuation method = ", method_impute)
    imputed = mice_missForest_impute(method_impute = method_impute, data = data, LABEL = LABEL, dirData = dirImpute, 
                                                 maxiter = maxiter, num_imputations = num_imputations, ntree = ntree,  
                                                visitOrder = visitOrder)
  
    MI_corr_W_means = vector("list",length = num_imputations)
    MI_corr_W_vars = vector("list",length = num_imputations)
    MI_corr = vector("list",length = num_imputations)

    for (s in 1: num_imputations){
      imp.num.n = imputed[[s]]
      LABEL = as.factor(imp.num.n[, "LABEL"])
      imp.num.n = imp.num.n[, -grep("LABEL", colnames(imp.num.n))]
      M_pearson = cor(imp.num.n, use = "pairwise.complete.obs", method = "pearson")
      M_kendall = cor(imp.num.n, use = "pairwise.complete.obs", method = "kendall")
      M_spearman = cor(imp.num.n, use = "pairwise.complete.obs", method = "spearman")
      MI_corr_W_means[[s]] = (M_pearson + M_kendall + M_spearman)/3   ## mean of s-th imputation 
      MI_corr_W_vars[[s]] = ((M_pearson-MI_corr_W_means[[s]])^2+
                            (M_kendall-MI_corr_W_means[[s]])^2+
                            (M_spearman-MI_corr_W_means[[s]])^2)/2    # between s-th imputation variance
      if (s== 1){
        MI_corr_global_mean = MI_corr_W_means[[s]]
        MI_corr_global_W_var =  MI_corr_W_vars[[s]]
      }else{
        MI_corr_global_mean = MI_corr_global_mean + MI_corr_W_means[[s]]
        MI_corr_global_W_var =  MI_corr_global_W_var + MI_corr_W_vars[[s]]
      }
       
    }
    
    MI_corr_global_mean = MI_corr_global_mean/s # global mean
    MI_corr_global_W_var = MI_corr_global_W_var/s   # global within-imputation correlation
    
    for (s in 1: num_imputations){
      if (s==1){
        MI_corr_B_var = (MI_corr_W_means[[s]] - MI_corr_global_mean)^2
      }else{
        MI_corr_B_var = MI_corr_B_var+(MI_corr_W_means[[s]] - MI_corr_global_mean)^2
      }
      
    }
    
    MI_corr_B_var = MI_corr_B_var/(s-1) # global between imputation variance
    MI_corr_global_var =  MI_corr_global_W_var + (1+1/s)*MI_corr_B_var
    #Wald test per le correlazioni
    wald_stats = abs(MI_corr_global_mean / MI_corr_global_var)
    p_vals_corrs = pchisq(wald_stats, df = 1,  lower.tail = FALSE)
    
    MI_corr_global_mean[!is.finite(p_vals_corrs)] = 0
    MI_corr_global_mean[p_vals_corrs>=thrConf] = 0
    corrplot::corrplot(MI_corr_global_mean, tl.cex = 0.65, tl.col = "black")
  }
  
#  correlazione con le LABEL
  
  for (method_impute in method_imputation){
    cat("start with impuation method = ", method_impute)
    imputed = mice_missForest_impute(method_impute = method_impute, data = data, LABEL = LABEL, dirData = dirImpute, 
                                     maxiter = maxiter, num_imputations = num_imputations, ntree = ntree,  
                                     visitOrder = visitOrder)
    
    MI_Labelcorrs = vector("list",length = num_imputations)
    
    for (s in 1: num_imputations){
      imp.num.n = imputed[[s]]
      LABEL = imp.num.n[, "LABEL"]
      imp.num.n = imp.num.n[, -grep("LABEL", colnames(imp.num.n))]
      MI_Labelcorrs[[s]] = rbind(apply(imp.num.n, 2,  cor, y = LABEL, use = "pairwise.complete.obs", method = "pearson"),
                        apply(imp.num.n, 2,  cor, y = LABEL, use = "pairwise.complete.obs", method = "kendall"),
                        apply(imp.num.n, 2,  cor, y = LABEL, use = "pairwise.complete.obs", method = "spearman"))
    }
  
    MI_Label_mean_corrs = rubin_rule_pool(MI_Labelcorrs)[["mean"]]
    MI_Label_var_corrs = rubin_rule_pool(MI_Labelcorrs)[["var"]]
    
    wald_stats = abs(MI_Label_mean_corrs / MI_Label_var_corrs)
    
    p_vals_corrs = pchisq(wald_stats, df = 1,  lower.tail = FALSE)
    
    MI_Label_mean_corrs [!is.finite(p_vals_corrs)] = 0
    MI_Label_mean_corrs [p_vals_corrs >= thrConf] = 0
    corrplot::corrplot(as.data.frame(MI_Label_mean_corrs), tl.cex = 0.65, tl.col = "black")
  }
  
 
 
}  #endfun




