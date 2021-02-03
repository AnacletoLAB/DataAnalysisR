main_clinical_data_analysis_NOSAT <- function(nIter = 5, ntree = 500, num_folds = 10, num.inner.folds = 7, num.boruta.folds = 5, strDir = '_covnet_score',  
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
  strdata = '_covnet_score_NOSAT'
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
  
  method_imputation = c("missForest")
  RF.results = vector("list", length = length(method_imputation))
  learner.results = vector("list", length = length(method_imputation))
  glm.results = vector("list", length = length(method_imputation))
  
 
  for (method_impute in method_imputation){
    cat("start with impuation method = ", method_impute)
    imputed = mice_missForest_impute(method_impute = method_impute, data = data, LABEL = LABEL, dirData = dirImpute, 
                                                 maxiter = maxiter, num_imputations, ntree = ntree,  
                                                visitOrder = visitOrder)
  
    cm.RF = vector("list",length = num_imputations)
    importances = vector("list", length = num_imputations)
    cm.learner =  vector("list",length = num_imputations)
    cm.glmnet =  vector("list",length = num_imputations)
    coeffs.glmnet = vector("list", length = num_imputations)
    
    for (s in 1: num_imputations){
      imp.num.n = imputed[[s]]
      LABEL = as.factor(imp.num.n[, "LABEL"])
      imp.num.n = imp.num.n[, -grep("LABEL", colnames(imp.num.n))]
      imp.n = transform.data.with.type(as.data.frame(imp.num.n))
      
      for (n.it in 1:nIter){
        ext.folds = folds.list[[n.it]]
        
        
        predictions.all <- RF_on_RF_and_Boruta(data = imp.n,LABEL = LABEL, ext.folds = ext.folds,  ntree = ntree,
                                               num.boruta.folds = num.boruta.folds, 
                                               num.inner.folds = num.inner.folds, num_folds = num_folds)
        
        
        # preds.glmnet <- glmnet_regression(data = imp.num.n, LABEL =LABEL, ext.folds, n_alpha = 31)
        
    
        if (n.it ==1){
          LABELS = predictions.all[["predictions.RF"]][,"LABELS"]
          predictions.RF = predictions.all[["predictions.RF"]][,"preds"]
          predictions.learner = predictions.all[["predictions.learner"]][,"preds"]
          probs.RF =  predictions.all[["predictions.RF"]][,"probs"]
          importances[[s]] = predictions.all[["RF.importances"]]
          
          # LABELS.glmnet =  preds.glmnet[["LABELS"]]
          # predictions.glmnet =  preds.glmnet[["predictions"]]
          # probs.glmnet =  preds.glmnet[["probs"]]
          # coeffs.glmnet[[s]] =  preds.glmnet[["coeffs"]]
          
          cm.RF[[s]] <- c (multiclass.confusion.matrix(as.factor(predictions.all[["predictions.RF"]][,"preds"]), 
                                                        as.factor(predictions.all[["predictions.RF"]][,"LABELS"]), metrics)[[2]][[4]],
                     "mcc" = MCC(as.factor(predictions.all[["predictions.RF"]][,"preds"]), 
                                       as.factor(predictions.all[["predictions.RF"]][,"LABELS"])), 
                     "auc" = pROC::auc(as.factor(predictions.all[["predictions.RF"]][,"LABELS"]), 
                                predictions.all[["predictions.RF"]]["probs"]$probs))
          
          cm.learner[[s]] <-  c(multiclass.confusion.matrix(as.factor(predictions.all[["predictions.learner"]][,"preds"]), 
                                                             as.factor(predictions.all[["predictions.learner"]][,"LABELS"]), metrics)[[2]][[4]],
                         "mcc"=  MCC(as.factor(predictions.all[["predictions.learner"]][,"preds"]), 
                                            as.factor(predictions.all[["predictions.learner"]][,"LABELS"])),
                         "auc" = pROC::auc( as.factor(predictions.all[["predictions.learner"]][,"LABELS"]),
                                      predictions.all[["predictions.learner"]]["preds"]$preds ))
                         
          
          # cm.glmnet[[s]] <- c(multiclass.confusion.matrix(as.factor(preds.glmnet[["LABELS"]]), 
          #                                                   as.factor(preds.glmnet[["predictions"]]), metrics)[[2]][[4]],
          #                   "mcc"= MCC(as.factor(preds.glmnet[["LABELS"]]), 
          #                                  as.factor(preds.glmnet[["predictions"]])),
          #                   "auc" = pROC::auc(as.factor(preds.glmnet[["LABELS"]]), preds.glmnet[["probs"]]))
          # 
          # 
        }else{
          LABELS = cbind(LABELS, predictions.all[["predictions.RF"]][,"LABELS"])
          predictions.RF = cbind(predictions.RF, predictions.all[["predictions.RF"]][,"preds"])
          predictions.learner = cbind(predictions.learner, predictions.all[["predictions.learner"]][,"preds"])
          probs.RF =  cbind(probs.RF, predictions.all[["predictions.RF"]][,"probs"])
          importances[[s]] = rbind(importances[[s]], predictions.all[["RF.importances"]])
          
          
          # LABELS.glmnet =  cbind(LABELS.glmnet, preds.glmnet[["LABELS"]])
          # predictions.glmnet =  cbind(predictions.glmnet, preds.glmnet[["predictions"]])
          # probs.glmnet =  cbind(probs.glmnet, preds.glmnet[["probs"]])
          # coeffs.glmnet[[s]] =  rbind(coeffs.glmnet[[s]], preds.glmnet[["coeffs"]])
          
          cm.RF[[s]] <- rbind(cm.RF[[s]] , c (multiclass.confusion.matrix(as.factor(predictions.all[["predictions.RF"]][,"preds"]), 
                                                                          as.factor(predictions.all[["predictions.RF"]][,"LABELS"]), metrics)[[2]][[4]],
                                              "mcc" = MCC(as.factor(predictions.all[["predictions.RF"]][,"preds"]), 
                                                          as.factor(predictions.all[["predictions.RF"]][,"LABELS"])), 
                                              "auc" = pROC::auc(as.factor(predictions.all[["predictions.RF"]][,"LABELS"]), 
                                                          predictions.all[["predictions.RF"]]["probs"]$probs)))
          
          
          cm.learner[[s]] <- rbind(cm.learner[[s]], c(multiclass.confusion.matrix(as.factor(predictions.all[["predictions.learner"]][,"preds"]), 
                                                                                  as.factor(predictions.all[["predictions.learner"]][,"LABELS"]), metrics)[[2]][[4]],
                                                      "mcc"=  MCC(as.factor(predictions.all[["predictions.learner"]][,"preds"]), 
                                                                  as.factor(predictions.all[["predictions.learner"]][,"LABELS"])),
                                                      "auc" = pROC::auc( as.factor(predictions.all[["predictions.learner"]][,"LABELS"]),
                                                                   predictions.all[["predictions.learner"]]["preds"]$preds )))
          
          # cm.glmnet[[s]] <- rbind(cm.glmnet[[s]] ,c(multiclass.confusion.matrix(as.factor(preds.glmnet[["LABELS"]]), 
          #                                                                       as.factor(preds.glmnet[["predictions"]]), metrics)[[2]][[4]],
          #                                           "mcc"= MCC(as.factor(preds.glmnet[["LABELS"]]), 
          #                                                      as.factor(preds.glmnet[["predictions"]])),
          #                                           "auc" = pROC::auc(as.factor(preds.glmnet[["LABELS"]]), preds.glmnet[["probs"]])))
          
        }
        
        
        #clear al, the variables use for this round
        
      
        predictions.all = NULL
     #   preds.glmnet = NULL
      }
    }
    
    RF.results[[grep(method_impute, method_imputation)]] = list("LABELS.RF" = LABELS, 
                                                        "cm.RF" = cm.RF, 
                                                        "predictions.RF" = predictions.RF, 
                                                        "probs.RF" = probs.RF, "importances.RF" = importances
                                                       )
    learner.results[[grep(method_impute, method_imputation)]] = list("LABELS.learner" = LABELS, 
                                                                "cm.learner" = cm.learner,
                                                                "predictions.learner" = predictions.learner
                                                                )
    # 
    # glm.results[[grep(method_impute, method_imputation)]] = list( "LABELS.glmnet" = LABELS.glmnet,
    #                                                                  "cm.glmnet" = cm.glmnet,
    #                                                                  "predictions.glmnet" = predictions.glmnet,
    #                                                                  "probs.glmnet" = probs.glmnet,
    #                                                                  "coeffs.glmnet" = coeffs.glmnet)
    
                                                        
    
    nameDf = file.path(nameResults, paste("results_", method_impute, "_", visitOrder, ".Rda", sep =""))
    save.image(file = nameDf)
    
  }    
  
  
 

  # feature importance plots
  for (nmethod in 1: length(method_imputation)){
    RF_results_mean_cv = vector("list", length = num_imputations)
#    Glm_results_mean_cv = vector("list", length = num_imputations)
    
   
    # normalizzo in modo che tutte le righe sommi a uno (ogni riga contiene le performance di tutte le variabili in un giro di fold cross validation)
    for (nlist in 1:length(RF.results[[nmethod]][["importances.RF"]])){
      
      for(nr in 1:nrow(RF.results[[nmethod]][["importances.RF"]][[nlist]])){ 
        idx = which(RF.results[[nmethod]][["importances.RF"]][[nlist]][nr, ]<0)
        if (length(idx)>0){
          ln = RF.results[[nmethod]][["importances.RF"]][[nlist]][nr,]
          ln[idx] = 0
          RF.results[[nmethod]][["importances.RF"]][[nlist]][nr,] = ln
        }
        RF.results[[nmethod]][["importances.RF"]][[nlist]][nr, ] = 
          RF.results[[nmethod]][["importances.RF"]][[nlist]][nr, ]/ sum(abs(RF.results[[nmethod]][["importances.RF"]][[nlist]][nr, ]))
        
      }
      RF_results_mean_cv[[nlist]] = as.data.frame(t(apply(RF.results[[nmethod]][["importances.RF"]][[nlist]][1:num_folds, ], 2, mean)))
     # Glm_results_mean_cv[[nlist]] = as.data.frame(t(apply(glm.results[[nmethod]][["coeffs.glmnet"]][[nlist]][1:num_folds, ], 2, mean)))
      for (nR in seq(from = 2, to = as.integer(nrow(RF.results[[nmethod]][["importances.RF"]][[nlist]])/num_folds))){
        ridx = seq(from = ((nR-1)*num_folds)+1 , to = nR*num_folds)
        #cat("mean over rows ", ridx, '\n')
        RF_results_mean_cv[[nlist]] = rbind(RF_results_mean_cv[[nlist]], 
                                            t(apply(RF.results[[nmethod]][["importances.RF"]][[nlist]][ridx, ], 2, mean) ))
      #  Glm_results_mean_cv[[nlist]] = rbind(Glm_results_mean_cv[[nlist]], 
       #                                      t(apply(glm.results[[nmethod]][["coeffs.glmnet"]][[nlist]][ridx, ], 2, mean) ))
      }
        
    } 
    
    RF_imp_mean = as.data.frame(rubin_rule_pool(RF_results_mean_cv)[["mean"]])
    colnames(RF_imp_mean) = "mean.importance"
    
    RF_imp_var = as.data.frame(rubin_rule_pool(RF_results_mean_cv)[["var"]])
    colnames(RF_imp_var) = "variance"
    
    RF_imp_se = as.data.frame(rubin_rule_pool(RF_results_mean_cv)[["se"]])
    colnames(RF_imp_se) = "stderr"
    
    # glm_coeff_mean = as.data.frame(rubin_rule_pool(Glm_results_mean_cv)[["mean"]])
    # colnames(glm_coeff_mean) = "mean.coeff"
    # 
    # glm_coeff_var = as.data.frame(rubin_rule_pool(Glm_results_mean_cv)[["var"]])
    # colnames(glm_coeff_var) = "variance"
    # 
    # glm_coeff_se = as.data.frame(rubin_rule_pool(Glm_results_mean_cv)[["se"]])
    # colnames(glm_coeff_se) = "stderr"
    
   # pchisq(rubin_rule_pool(q, df,  lower.tail = FALSE)
    q_RF = RF_imp_mean/RF_imp_var
    colnames(q_RF) = "wald"
    p_val_RF = as.data.frame(pchisq(as.numeric(q_RF$wald), df = 1,  lower.tail = FALSE))
    colnames(p_val_RF ) = "pval"
    
    
    # q_glm = glm_coeff_mean/glm_coeff_var
    # colnames(q_glm) = "wald"
    # p_val_glm = as.data.frame(pchisq(as.numeric(q_glm$wald), df = 1,  lower.tail = FALSE))
    # colnames(p_val_glm ) = "pval"
    # 
    if (nmethod ==1){
      df_RF_importances = cbind("imputation" = method_imputation[nmethod], 
                                "variable"= factor(row.names( RF_imp_mean), levels = row.names( RF_imp_mean)),  
                                RF_imp_mean, RF_imp_var, RF_imp_se, p_val_RF)
      # df_glm_importances = cbind("imputation" = method_imputation[nmethod],
      #                            "variable"= factor(row.names( glm_coeff_mean), levels = row.names(  glm_coeff_mean)), 
      #                            glm_coeff_mean, glm_coeff_var, glm_coeff_se, p_val_glm)
    }else{
      df_RF_importances = rbind(df_RF_importances, cbind("imputation" = method_imputation[nmethod],
                                                         "variable"= factor(row.names( RF_imp_mean), levels = row.names( RF_imp_mean)), 
                                                         RF_imp_mean, RF_imp_var, RF_imp_se, p_val_RF))
      # df_glm_importances = rbind(df_glm_importances, cbind("imputation" = method_imputation[nmethod], 
      #                                                      "variable"= factor(row.names( glm_coeff_mean), levels = row.names(  glm_coeff_mean)), 
      #                                                      glm_coeff_mean, glm_coeff_var, glm_coeff_se, p_val_glm))
    }
  }
  
  
  idx_del_RF = which(is.na(df_RF_importances$pval) | (df_RF_importances$pval>thrConf))
  if (length(idx_del_RF)>0){
    df_RF_importances[idx_del_RF, c("mean.importance", "variance", "stderr")]=0
  }
    
  pRF<- ggplot(df_RF_importances, aes(x = variable,  y=mean.importance, fill=imputation)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
     geom_errorbar(aes(ymin=mean.importance-stderr, ymax=mean.importance+stderr), width=.2,
                   position=position_dodge(.9)) 
    # geom_errorbar(aes(ymin=mean.importance-sqrt(variance)/sqrt(num_imputations), ymax=mean.importance+sqrt(variance)/sqrt(num_imputations)), width=.2,
    #               position=position_dodge(.9)) 
    
  # Finished bar plot
  pRF = pRF+labs(title="RF Importance", x="Variables", y = "RF Importance")+
    theme(axis.text.x = element_text(angle = 90))
  print(pRF)
  
  
  # idx_del_Glm = which(is.na(df_glm_importances$pval) | (df_glm_importances$pval>thrConf))
  # if (length(idx_del_RF)>0){
  #   df_glm_importances[idx_del_Glm, c("mean.coeff", "variance", "stderr")]=0
  # }
  # pGLM<- ggplot(df_glm_importances, aes(x = variable,  y=mean.coeff, fill=imputation)) + 
  #   geom_bar(stat="identity", color="black", 
  #            position=position_dodge()) +
  #   geom_errorbar(aes(ymin=mean.coeff-stderr, ymax=mean.coeff+stderr), width=.2,
  #                  position=position_dodge(.9)) 
  #     # geom_errorbar(aes(ymin=mean.coeff-sqrt(variance)/sqrt(num_imputations), ymax=mean.coeff+sqrt(variance)/sqrt(num_imputations)), width=.2,
  #   #               position=position_dodge(.9)) 
  # 
  # pGLM = pGLM+labs(title="glm Importance", x="Variables", y = "GLM coefficient")+
  #   theme(axis.text.x = element_text(angle = 90))
  # print(pGLM)
  # 
  # grid.arrange(pRF, pGLM, nrow = 2)

  df_RF_missForest = df_RF_importances[which(df_RF_importances$imputation=="missForest"), ]
  # df_glm_missForest = df_glm_importances[which(df_glm_importances$imputation=="missForest"), ]
  # 
  mega_df = cbind("RF.imp"= df_RF_missForest$mean.importance, "RF.stderr" = df_RF_missForest$stderr,
                  "RF.pval"= df_RF_missForest$pval)
                  # "glm.coeff" =  df_glm_missForest$mean.coeff, "glm.stderr" =  df_glm_missForest$stderr, 
                  # "glm.pval" =  df_glm_missForest$pval) 
  row.names(mega_df) = levels(df_RF_missForest$variable)
  
  write.csv(mega_df, file = file.path(nameResults, "importances_RF_NOSAT.csv"))
    
  #--------------- END feature importance plots
  
  cm_all = vector("list", length(method_imputation))
  # wilcoxon test for comparing AUC sensitivity specificity, etc
  risk.models = factor(c("RF", "AT"), levels = c("RF", "AT"))
  cnames = colnames(RF.results[[nmethod]][["cm.RF"]][[1]])
  cnames[which(cnames == "Balanced Accuracy")] = "Accuracy"
  for (nmethod in 1: length(method_imputation)){
    for (nimp in 1: num_imputations){
      df = rbind(cbind("risk.model" = risk.models[1] , as.data.frame(t(apply(RF.results[[nmethod]][["cm.RF"]][[nimp]], 2, mean)))), 
                 cbind("risk.model" = risk.models[2] , as.data.frame(t(apply(learner.results[[nmethod]][["cm.learner"]][[nimp]], 2, mean) ))))
      colnames(df) = c("risk.model", cnames)
      if (nimp ==1){
        cm_all[[nmethod]] = df
      }else{
        cm_all[[nmethod]] = rbind(cm_all[[nmethod]], df)
      }
    }
  }
  results_all = rbind(cbind("imputation" = method_imputation[1],  cm_all[[1]]),
                   cbind("imputation" = method_imputation[2],  cm_all[[2]]))
  
  p_mice_missRF = rbind ( cbind( "less" = wilcox.test(auc ~ imputation, data = results_all, paired = FALSE, alternative = "less")[["p.value"]],
                                "greater" = wilcox.test(auc ~ imputation, data = results_all, paired = FALSE, alternative = "greater")[["p.value"]],
                                "less.RF" = wilcox.test(auc ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "less")[["p.value"]],
                                "greater.RF" = wilcox.test(auc ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                "less.AT" = wilcox.test(auc ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "less")[["p.value"]],
                                "greater.AT" = wilcox.test(auc ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                "less.Glm" = wilcox.test(auc ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "less")[["p.value"]],
                                "greater.Glm" = wilcox.test(auc ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "greater")[["p.value"]]),
                  cbind("less" = wilcox.test(Sensitivity ~ imputation, data = results_all, paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater" = wilcox.test(Sensitivity ~ imputation, data = results_all, paired = FALSE, alternative = "greater")[["p.value"]],
                                       "less.RF" = wilcox.test(Sensitivity ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater.RF" = wilcox.test(Sensitivity ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                       "less.AT" = wilcox.test(Sensitivity ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater" = wilcox.test(Sensitivity ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                       "less" = wilcox.test(Sensitivity ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater" = wilcox.test(Sensitivity ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "greater")[["p.value"]]),
                 cbind("less" = wilcox.test(Specificity ~ imputation, data = results_all, paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater" = wilcox.test(Specificity ~ imputation, data = results_all, paired = FALSE, alternative = "greater")[["p.value"]],
                                       "less.RF" = wilcox.test(Specificity ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater.RF" = wilcox.test(Specificity ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                       "less.AT" = wilcox.test(Specificity ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater" = wilcox.test(Specificity ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                       "less" = wilcox.test(Specificity ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "less")[["p.value"]],
                                       "greater" = wilcox.test(Specificity ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "greater")[["p.value"]]),
                 cbind("less" = wilcox.test(F1 ~ imputation, data = results_all, paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater" = wilcox.test(F1 ~ imputation, data = results_all, paired = FALSE, alternative = "greater")[["p.value"]],
                                    "less.RF" = wilcox.test(F1 ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater.RF" = wilcox.test(F1 ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                    "less.AT" = wilcox.test(F1 ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater" = wilcox.test(F1 ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                    "less" = wilcox.test(F1 ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater" = wilcox.test(F1 ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "greater")[["p.value"]]),
                  cbind("less" = wilcox.test(Accuracy ~ imputation, data = results_all, paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater" = wilcox.test(Accuracy ~ imputation, data = results_all, paired = FALSE, alternative = "greater")[["p.value"]],
                                    "less.RF" = wilcox.test(Accuracy ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater.RF" = wilcox.test(Accuracy ~ imputation, data = results_all[which(results_all$risk.model=="RF"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                    "less.AT" = wilcox.test(Accuracy ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater" = wilcox.test(Accuracy ~ imputation, data = results_all[which(results_all$risk.model=="AT"),], paired = FALSE, alternative = "greater")[["p.value"]],
                                    "less" = wilcox.test(Accuracy ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "less")[["p.value"]],
                                    "greater" = wilcox.test(Accuracy ~ imputation, data = results_all[which(results_all$risk.model=="Glm"),], paired = FALSE, alternative = "greater")[["p.value"]]))
                 
  rownames(p_mice_missRF) = c("AUC", "Sensitivity", "Specificity", "F1", "Accuracy" )
  
  write.csv(p_mice_missRF, file = file.path(nameResults, 'wilcox.compare_miceRF_missForest.csv'))
  
  # compare risk predictors
  results_missForest = results_all[which(results_all$imputation=="missForest"), ]
  results_missForest_RF_Glm = results_missForest[-which(results_missForest$risk.model == "AT"), ]
  results_missForest_RF_AT = results_missForest[-which(results_missForest$risk.model == "Glm"), ]
  
  missForest_RF_AT_GLM = rbind(cbind("RF vs AT" = wilcox.test(auc ~ risk.model, data = results_missForest_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
        "RF vs Glm" = wilcox.test(auc ~ risk.model, data = results_missForest_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
  cbind("RF vs AT" = wilcox.test(Sensitivity ~ risk.model, data = results_missForest_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
        "RF vs Glm" = wilcox.test(Sensitivity ~ risk.model, data = results_missForest_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
  cbind("RF vs AT" = wilcox.test(Specificity ~ risk.model, data = results_missForest_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
        "RF vs Glm" = wilcox.test(Specificity ~ risk.model, data = results_missForest_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
  cbind("RF vs AT" = wilcox.test(F1 ~ risk.model, data = results_missForest_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
        "RF vs Glm" = wilcox.test(F1 ~ risk.model, data = results_missForest_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
  cbind("RF vs AT" = wilcox.test(Accuracy ~ risk.model, data = results_missForest_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
        "RF vs Glm" = wilcox.test(Accuracy ~ risk.model, data = results_missForest_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]))

  rownames( missForest_RF_AT_GLM )  =  c("AUC", "Sensitivity", "Specificity", "F1", "Accuracy" )
  
  results_miceRF = results_all[which(results_all$imputation=="miceRF"), ]
  results_miceRF_RF_Glm = results_missForest[-which(results_miceRF$risk.model == "AT"), ]
  results_miceRF_RF_AT = results_missForest[-which(results_miceRF$risk.model == "Glm"), ]
  
  miceRF_RF_AT_GLM = rbind(cbind("RF vs AT" = wilcox.test(auc ~ risk.model, data = results_miceRF_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
                                     "RF vs Glm" = wilcox.test(auc ~ risk.model, data = results_miceRF_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
                               cbind("RF vs AT" = wilcox.test(Sensitivity ~ risk.model, data = results_miceRF_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
                                     "RF vs Glm" = wilcox.test(Sensitivity ~ risk.model, data = results_miceRF_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
                               cbind("RF vs AT" = wilcox.test(Specificity ~ risk.model, data = results_miceRF_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
                                     "RF vs Glm" = wilcox.test(Specificity ~ risk.model, data = results_miceRF_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
                               cbind("RF vs AT" = wilcox.test(F1 ~ risk.model, data = results_miceRF_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
                                     "RF vs Glm" = wilcox.test(F1 ~ risk.model, data = results_miceRF_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]),
                               cbind("RF vs AT" = wilcox.test(Accuracy ~ risk.model, data = results_miceRF_RF_AT, paired = FALSE, alternative = "greater")[["p.value"]], 
                                     "RF vs Glm" = wilcox.test(Accuracy ~ risk.model, data = results_miceRF_RF_Glm, paired = FALSE, alternative = "greater")[["p.value"]]))
  
  rownames(miceRF_RF_AT_GLM )  =  c("AUC", "Sensitivity", "Specificity", "F1", "Accuracy" )
  
  comp_RF_GLM_AT = rbind(cbind("imputation" = "missForest", t( missForest_RF_AT_GLM)), cbind("imputation" = "miceRF", t(miceRF_RF_AT_GLM)))
  
  write.csv(comp_RF_GLM_AT, file = file.path(nameResults, 'wilcox.compare_risk_predictors.csv'))
  
  
  # rubin's rules for AUC sensitivity specificity, etc
  
  RF_imp_auc = matrix(0,ncol = length(method_imputation), nrow = 4)
  learner_imp_auc = matrix(0,ncol = length(method_imputation), nrow = 4)
  gml_imp_auc = matrix(0,ncol = length(method_imputation), nrow = 4)
  
  RF_imp_sensitivity = matrix(0,ncol = length(method_imputation), nrow = 4)
  learner_imp_sensitivity = matrix(0,ncol = length(method_imputation), nrow = 4)
  gml_imp_sensitivity = matrix(0,ncol = length(method_imputation), nrow = 4)
  
  RF_imp_specificity = matrix(0,ncol = length(method_imputation), nrow = 4)
  learner_imp_specificity = matrix(0,ncol = length(method_imputation), nrow = 4)
  gml_imp_specificity = matrix(0,ncol = length(method_imputation), nrow = 4)
  
  RF_imp_accuracy = matrix(0,ncol = length(method_imputation), nrow = 4)
  learner_imp_accuracy = matrix(0,ncol = length(method_imputation), nrow = 4)
  gml_imp_accuracy = matrix(0,ncol = length(method_imputation), nrow = 4)
  
  RF_imp_F1 = matrix(0,ncol = length(method_imputation), nrow = 4)
  learner_imp_F1 = matrix(0,ncol = length(method_imputation), nrow = 4)
  gml_imp_F1 = matrix(0,ncol = length(method_imputation), nrow = 4)
  
  prec_round = 5

  
  for (nmethod in 1: length(method_imputation)){
    RF_imp_auc[1, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["mean"]]["auc"],2)
    RF_imp_auc[2, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["var"]]["auc"],prec_round)
    RF_imp_auc[3, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["within"]]["auc"],prec_round)
    RF_imp_auc[4, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["between"]]["auc"],prec_round)
    
    RF_imp_sensitivity[1, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["mean"]]["Sensitivity"],2)
    RF_imp_sensitivity[2, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["var"]]["Sensitivity"],prec_round)
    RF_imp_sensitivity[3, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["within"]]["Sensitivity"],prec_round)
    RF_imp_sensitivity[4, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["between"]]["Sensitivity"],prec_round)
    
    RF_imp_specificity[1, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["mean"]]["Specificity"],2)
    RF_imp_specificity[2, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["var"]]["Specificity"],prec_round)
    RF_imp_specificity[3, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["within"]]["Specificity"],prec_round)
    RF_imp_specificity[4, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["between"]]["Specificity"],prec_round)
    
    RF_imp_accuracy[1, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["mean"]]["Balanced Accuracy"],2)
    RF_imp_accuracy[2, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["var"]]["Balanced Accuracy"],prec_round)
    RF_imp_accuracy[3, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["within"]]["Balanced Accuracy"],prec_round)
    RF_imp_accuracy[4, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["between"]]["Balanced Accuracy"],prec_round)
    
    RF_imp_F1[1, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["mean"]]["F1"],2)
    RF_imp_F1[2, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["var"]]["F1"],prec_round)
    RF_imp_F1[3, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["within"]]["F1"],prec_round)
    RF_imp_F1[4, nmethod] = round(rubin_rule_pool(RF.results[[nmethod]][["cm.RF"]])[["between"]]["F1"],prec_round)
    
    
    
  
    
    
    
    learner_imp_auc[1, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["mean"]]["auc"],2)
    learner_imp_auc[2, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["var"]]["auc"],prec_round)
    learner_imp_auc[3, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["within"]]["auc"],prec_round)
    learner_imp_auc[4, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["between"]]["auc"],prec_round)
    
    
    learner_imp_sensitivity[1, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["mean"]]["Sensitivity"],2)
    learner_imp_sensitivity[2, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["var"]]["Sensitivity"],prec_round)
    learner_imp_sensitivity[3, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["within"]]["Sensitivity"],prec_round)
    learner_imp_sensitivity[4, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["between"]]["Sensitivity"],prec_round)
    
    learner_imp_specificity[1, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["mean"]]["Specificity"],2)
    learner_imp_specificity[2, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["var"]]["Specificity"],prec_round)
    learner_imp_specificity[3, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["within"]]["Specificity"],prec_round)
    learner_imp_specificity[4, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["between"]]["Specificity"],prec_round)
    
    learner_imp_accuracy[1, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["mean"]]["Balanced Accuracy"],2)
    learner_imp_accuracy[2, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["var"]]["Balanced Accuracy"],prec_round)
    learner_imp_accuracy[3, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["within"]]["Balanced Accuracy"],prec_round)
    learner_imp_accuracy[4, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["between"]]["Balanced Accuracy"],prec_round)
    
    learner_imp_F1[1, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["mean"]]["F1"],2)
    learner_imp_F1[2, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["var"]]["F1"],prec_round)
    learner_imp_F1[3, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["within"]]["F1"],prec_round)
    learner_imp_F1[4, nmethod] = round(rubin_rule_pool(learner.results[[nmethod]][["cm.learner"]])[["between"]]["F1"],prec_round)
    
    
    
    gml_imp_auc[1, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["mean"]]["auc"],2)
    gml_imp_auc[2, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["var"]]["auc"],prec_round)
    gml_imp_auc[3, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["within"]]["auc"],prec_round)
    gml_imp_auc[4, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["between"]]["auc"],prec_round)
    
    gml_imp_sensitivity[1, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["mean"]]["Sensitivity"],2)
    gml_imp_sensitivity[2, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["var"]]["Sensitivity"],prec_round)
    gml_imp_sensitivity[3, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["within"]]["Sensitivity"],prec_round)
    gml_imp_sensitivity[4, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["between"]]["Sensitivity"],prec_round)
    
    gml_imp_specificity[1, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["mean"]]["Specificity"],2)
    gml_imp_specificity[2, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["var"]]["Specificity"],prec_round)
    gml_imp_specificity[3, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["within"]]["Specificity"],prec_round)
    gml_imp_specificity[4, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["between"]]["Specificity"],prec_round)
    
    gml_imp_accuracy[1, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["mean"]]["Balanced Accuracy"],2)
    gml_imp_accuracy[2, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["var"]]["Balanced Accuracy"],prec_round)
    gml_imp_accuracy[3, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["within"]]["Balanced Accuracy"],prec_round)
    gml_imp_accuracy[4, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["between"]]["Balanced Accuracy"],prec_round)
    
    gml_imp_F1[1, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["mean"]]["F1"],2)
    gml_imp_F1[2, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["var"]]["F1"],prec_round)
    gml_imp_F1[3, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["within"]]["F1"],prec_round)
    gml_imp_F1[4, nmethod] = round(rubin_rule_pool(glm.results[[nmethod]][["cm.glmnet"]])[["between"]]["F1"],prec_round)
    
    important_vars = sort(rubin_rule_pool(glm.results[[nmethod]][["coeffs.glmnet"]])[["mean"]], decreasing = TRUE)
  }
  
  
  nameDf = file.path(nameResults, paste("results_all_", visitOrder, "_m", num_imputations, ".Rda", sep =""))
  save.image(file = nameDf)
  
  
  
 
 
}  #endfun




