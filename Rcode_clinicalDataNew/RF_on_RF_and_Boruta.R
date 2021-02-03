RF_on_RF_and_Boruta <- function(data,LABEL, ext.folds, ntree = 31, num_folds = 10, num.inner.folds = 5, num.boruta.folds = 3) {

    
  library(randomForest)
  library(randomForestExplainer)
  library(inTrees)
  library(caret)
  
  ntree.rules.extract=ntree
  fraction.node.size = 30
  depth.rules.to.extract = 11
  source(file.path('.', 'RF_on_Boruta_new.R'))
  source(file.path('.', 'tune_xgb.R'))
  source(file.path('.', 'boruta_n_fold_select.R'))
  
 # ext.folds <- createFolds(y= LABEL, k = num_folds, list = TRUE, returnTrain = FALSE)

  
  imp.rf = data.frame(matrix(0, nrow = length(ext.folds), ncol = ncol(data)))
  colnames(imp.rf) = colnames(data)
  
  num_classes = length(levels(LABEL))
  for (num_f in 1:length(ext.folds)){
    cat('num external fold = ', num_f, '\n')
    ns = ext.folds[[num_f]]
    # # apro il file dei tipi in cui le features non sono ancora state selezionate da Boruta: cosÃ¬ ho i nomi di tutte le features
    # # butto via le colonne TAG che non corrispondono a features
    data.train = data[-ns,]
    lab.train = LABEL[-ns]
    
    data.test = data[ns,]
    lab.test = LABEL[ns]
    
    cat("In RF RF_on_RF_and_Boruta LABEL is categorical before RF_and_Boruta? ", is.factor(lab.train), '\n')
    
    result <- RF_on_Boruta_new(data = data.train, LABEL = lab.train, ntree=ntree, num_folds = num.inner.folds, num.boruta.folds = num.boruta.folds) 
    
    importances = result[["importances"]]
    preds = result[["predictions"]]
    labs = result[["LABELS"]]
    
  
    var_important = as.data.frame(apply(importances, 1, mean))
    colnames(var_important) = "mean.Imp"
    
    var_important = as.data.frame(apply(var_important, 2, sort, decreasing=TRUE))
   
    sorted = cbind(var_important,  "cumsum" = round(cumsum(var_important$mean.Imp/sum(var_important$mean.Imp)),2))
    sel_vars = unlist(rownames(sorted[which(sorted$cumsum<=0.95), ]))
    del_vars = unlist(rownames(sorted[which(sorted$cumsum>0.95), ]))
    
    if (length(sel_vars)==0){
      warning("NO VARIABLES SELECTED IN EXTERNAL FOLD!")
      sel_vars = colnames(data.train)
      del_vars = NULL
    }
    cat("\n\n", length(sel_vars), " variables selected in external fold  ",num_f,  "\n")
    cat(sel_vars)
    cat("\n")
    
    # non mi serve la LABEL nel data frame perchè poi la toglielvo sempre!
    data.train = data.train[,sel_vars]
    cat(colnames(data.train))
    
    data.test =  data.test[,sel_vars]
    cat(colnames(data.test))
    
    cat("\n", "selected vars = ", colnames(data.train))
    cat("\n")
    
    #tune_xgb(data.train, nameLab)
# 
#     trC <- trainControl(method = "repeatedcv", search= "random")
#     
#     
#     glm_fit <- caret::train(LABEL ~ ., data = data.train,
#                          method = "glmboost",
#                          family = "binomial",
#                          ## Create 20 random parameter values
#                          tuneLength = 20,
#                          trControl = trC)
#     
#     
#     
#     
#     preds =  predict(glm_fit, data.test[,-grep(nameLab, colnames(data.test))])
#     if (num_f ==1){
#       predictions.glm = as.data.frame(cbind("preds" = as.factor(preds), "LABELS" = as.factor(lab.test))) 
#     }else{
#       predictions.glm = rbind(predictions.glm, 
#                               as.data.frame(cbind("preds" = as.factor(preds), "LABELS" = as.factor(lab.test))) )
#     }
    
    # nnetFit <- train(LABEL ~ ., data = data.train,
    #                  method = "nnet",
    #                  preProcess = "range",
    #                  tuneLength = 2,
    #                  trace = FALSE,
    #                  maxit = 100)
    # coef(nnetFit$finalModel)
    # preds =  predict(nnetFit, data.test[,-grep(nameLab, colnames(data.test))])
    # if (num_f ==1){
    #   predictions.net = as.data.frame(cbind("preds" = as.factor(preds), "LABELS" = as.factor(lab.test))) 
    # }else{
    #   predictions.net = rbind(predictions.net, 
    #                       as.data.frame(cbind("preds" = as.factor(preds), "LABELS" = as.factor(lab.test))) )
    # }
    # 
    
    clw <- rep(1 / length(lab.train),num_classes)
    names(clw) = levels(lab.train)
    spsz <- rep(min(table(lab.train)),num_classes)
    names(spsz) = levels(lab.train)
    
    for (nc in levels(lab.train)){ clw[nc] = clw[nc]*length(which(lab.train == nc))}
    
    rf_classifier = tuneRF(y = lab.train, x=data.train,importance=TRUE, classwt = clw,
                           strata = lab.train, sampsize = spsz, nodesize = nrow(data.train)/ fraction.node.size,
                           mtryStart = ncol(data.train)/2, ntreeTry=ntree, stepFactor=1, improve=0.01,trace=TRUE, 
                           plot=FALSE, doBest=TRUE)
    
    cat("rf is for classification? ", rf_classifier$type, '\n')
    
    imps = importance( rf_classifier, type = 1)
    if (any(imps<0)){
      imps[which(imps<0)] = 0
    }
    imp.rf[num_f, sel_vars] = imps 
    
    preds.prob <-predict(rf_classifier,data.test, type = "prob")
    preds <- as.factor(predict(rf_classifier,data.test))
    
    if (num_f ==1){
      predictions = as.data.frame(cbind("preds" = as.factor(preds), "probs" = preds.prob[,"1"], "LABELS" = as.factor(lab.test))) 
    }else{
      predictions = rbind(predictions, 
                          as.data.frame(cbind("preds" = as.factor(preds),"probs" = preds.prob[,"1"], "LABELS" = as.factor(lab.test))) )
    }
    
    treeList <- inTrees::RF2List(rf_classifier)
    ruleExec <- unique(inTrees::extractRules(treeList,data.train, digits=2, ntree = ntree, maxdepth = depth.rules.to.extract))
    ruleMetric <- inTrees::getRuleMetric(ruleExec,data.train,lab.train) # measure rules
    ruleMetric <- inTrees::pruneRule(ruleMetric,data.train,lab.train) # prune each rule
    
    # create the learner for this fold
    learner <- inTrees::buildLearner(ruleMetric,data.train,lab.train)
    preds <- inTrees::applyLearner(learner,data.test)
    if (num_f ==1){
      predictions.learner = as.data.frame(cbind("preds" = as.factor(preds), "LABELS" = as.factor(lab.test))) 
    }else{
      predictions.learner = rbind(predictions.learner, 
                          as.data.frame(cbind("preds" = as.factor(preds), "LABELS" = as.factor(lab.test))) )
    }
    
    
    ruleExec <- as.matrix(ruleMetric[, "condition"])
     # since features are different for each fold, they must be converted to be merged at the end
    ruleExec <- convert.rules(ruleExec = ruleExec, new.idx = match(colnames(data.train), colnames(data)))
    ruleMetric[, "condition"] = ruleExec
    if (num_f==1){   
     all.rules = unique(ruleExec)
     all.ruleMetric = unique(ruleMetric)
    }else{
      all.rules = unique(rbind(all.rules, ruleExec))
      all.ruleMetric = unique(rbind(all.ruleMetric, ruleMetric))
    }
  }  
   
  return(list("RF.importances" = imp.rf, "predictions.RF" = predictions, "predictions.learner" = predictions.learner, 
              "rules" = all.rules, "rules.Metric" = all.ruleMetric))
 
}  #endfun