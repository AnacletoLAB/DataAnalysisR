boruta_n_fold_select <- function(data.boruta, label.boruta, max.boruta.runs = 101, num_folds = 10, select.tentative = TRUE){
# creting num_folds folds to insert randomness
  # on each repetition take 9 folds and apply Boruta
  # keep also 'tentative' variables if select.tentative == TRUE
  #  keep only 'confirmed' variables if select.tentative == FALSE
  
  library(Boruta)

  invotes = rep(0,times = ncol(data.boruta))
  label.boruta = as.factor(label.boruta)
  if (num_folds ==1){
    BorutaOnDataTrain <- Boruta(data.boruta, label.boruta, maxRuns = max.boruta.runs)
    idx_sel = grep('Confirmed', BorutaOnDataTrain$finalDecision)
    if (select.tentative){ idx_sel <- c(idx_sel, grep('Tentative', BorutaOnDataTrain$finalDecision))}
    invotes[idx_sel] <- 1
  }else{
    inner.folds <- createFolds(label.boruta, k = num_folds, list = TRUE, returnTrain = FALSE)
    
    for (ff in inner.folds){
        cat('b ')
  
        #APPLY BORUTA
        BorutaOnDataTrain <- Boruta(data.boruta[-ff,], label.boruta[-ff], maxRuns = max.boruta.runs)
        idx_sel = grep('Confirmed', BorutaOnDataTrain$finalDecision)
        if (select.tentative){ idx_sel <- c(idx_sel, grep('Tentative', BorutaOnDataTrain$finalDecision))}
        invotes[idx_sel] <- invotes[idx_sel]+1
  
    }
    cat('\n')
    
  }
  
  return(invotes)

}
