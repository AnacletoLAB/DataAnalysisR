tune_xgb <- function(data, nameLab = "LABEL"){
  library(xgboost)

  cv_folds <- KFold(as.matrix(data[,grep(nameLab, colnames(data))]), nfolds = 5, 
                    stratified = TRUE, seed = 50)
  xgb.cv.bayes <- function(max.depth, min_child_weight, subsample, colsample_bytree, gamma){
    cv <- xgb.cv(params = list(booster = 'gbtree', eta = 0.05,
                               max_depth = max.depth,
                               min_child_weight = min_child_weight,
                               subsample = subsample,
                               colsample_bytree = colsample_bytree,
                               gamma = gamma,
                               lambda = 1, alpha = 0,
                               objective = 'binary:logistic',
                               metrics = 'aucpr'),
                 data = data.matrix(data[,-grep(nameLab, colnames(data))]),
                 label = as.matrix(data[, grep(nameLab, colnames(data))]),
                 nround = 500, folds = cv_folds, prediction = TRUE,
                 showsd = TRUE, early_stopping_rounds = 5, maximize = TRUE,
                 verbose = 0
    )
    list(Score = cv$dt[, max(aucpr)],
         Pred = cv$pred)
  }
  
  xgb.bayes.model <- BayesianOptimization(
    xgb.cv.bayes,
    bounds = list(max.depth = c(2L, 12L),
                  min_child_weight = c(1L, 10L),
                  subsample = c(0.5, 1),
                  colsample_bytree = c(0.1, 0.4),
                  gamma = c(0, 10)
    ),
    init_grid_dt = NULL,
    init_points = 10,  # number of random points to start search
    n_iter = 20, # number of iterations after initial random points are set
    acq = 'ucb', kappa = 2.576, eps = 0.0, verbose = TRUE
  )
  
  return(xgb.bayes.model)
  
}
  
