#
#df <- data.frame("Prediction" = predictions, "Reference" = ref.labels)

#da https://www.datascienceblog.net/post/machine-learning/performance-measures-multi-class-problems/


calculate.accuracy <- function(predictions, ref.labels) {
  return(length(which(predictions == ref.labels)) / length(ref.labels))


}


calculate.w.accuracy <- function(predictions, ref.labels, flag_print = TRUE) {
  lvls <- levels(ref.labels)
  
  
  accs <- lapply(lvls, function(x) {
    idx <- which(ref.labels == x)
    return(calculate.accuracy(predictions[idx], ref.labels[idx]))
  })
  acc <- mean(unlist(accs))
  if (flag_print){
      print(paste("weigthed accuracy is: ", round(acc, 2)))
  }
  return(acc)
}


#multiclass confusion matrix
#metrics contiene le metriche da mostrare
multiclass.confusion.matrix <- function(predictions, ref.labels, positive ="1", metrics= 'all', flag_print = TRUE){
    #metrics <- c("Precision", "Recall")
  
  library(caret) # for confusionMatrix function
  cm <- vector("list", length(levels(ref.labels)))
  for (i in seq_along(cm)) {
    positive.class <- levels(ref.labels)[i]
    # in the i-th iteration, use the i-th class as the positive class
    cm[[i]] <- confusionMatrix(as.factor(predictions), as.factor(ref.labels), 
                               positive = positive.class)
  }
  if (flag_print){
    if (length(metrics)==1){
      if (metrics == 'all'){
        print(cm[[1]]$byClass)
      }
    }else{
      if (length(levels(ref.labels))<=2){
        cm = cbind( cm[[1]]$byClass[metrics], "mcc" = MCC(predictions, ref.labels))
        if (flag_print){
          print(cm)
        }
      }else{ if (flag_print){cm = cm[[1]]$byClass[, metrics]}}
    }
  }
  
  return(cm)
}

MCC <- function(predictions, ref.labels){
  library(caret) # for confusionMatrix function
  num_classes = length(levels(ref.labels)) 
  for (i in 1:num_classes) {
    positive.class <- levels(ref.labels)[i]
    # in the i-th iteration, use the i-th class as the positive class
    conf.mat <- confusionMatrix(as.factor(predictions), as.factor(ref.labels), 
                               positive = positive.class)[[2]]
    TP = conf.mat[2,2]
    TN = conf.mat[1,1]
    FN = conf.mat[1,2]
    FP = conf.mat[2,1]
    if (i == 1){val_mcc = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    }else{ val_mcc= c(val_mcc, ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))}
  }
  if (num_classes ==2){ 
    cat("2 classes: ", levels(ref.labels)[1], " considered as positive class (mcc = ", val_mcc[2] , ")", "\n")
    return( val_mcc[2])
  }else{ return(unlist(val_mcc))}
  
}

get.conf.stats <- function(cm) {
  out <- vector("list", length(cm))
  for (i in seq_along(cm)) {
    x <- cm[[i]]
    tp <- x$table[x$positive, x$positive] 
    fp <- sum(x$table[x$positive, colnames(x$table) != x$positive])
    fn <- sum(x$table[colnames(x$table) != x$positie, x$positive])
    # TNs are not well-defined for one-vs-all approach
    elem <- c(tp = tp, fp = fp, fn = fn)
    out[[i]] <- elem
  }
  df <- do.call(rbind, out)
  rownames(df) <- unlist(lapply(cm, function(x) x$positive))
  return(as.data.frame(df))
}



#micro_averaged F1: prende in input la multiclass confusion matrix
get.micro.f1 <- function(cm, flag_print = TRUE) {
  cm.summary <- get.conf.stats(cm)
  tp <- sum(cm.summary$tp)
  fn <- sum(cm.summary$fn)
  fp <- sum(cm.summary$fp)
  pr <- tp / (tp + fp)
  re <- tp / (tp + fn)
  f1 <- 2 * ((pr * re) / (pr + re))
  if (flag_print){print(paste("Micro F1 is: ", round(f1, 2)))}
  return(f1)

  
}

#macro_averaged F1: prende in input la multiclass confusion matrix
get.macro.f1 <- function(cm, flag_print = TRUE) {
  c <- cm[[1]]$byClass # a single matrix is sufficient
  if (is.matrix(c)){
    re <- sum(c[, "Recall"]) / nrow(c)
    pr <- sum(c[, "Precision"]) / nrow(c)
  }else{
    re <- c["Recall"]
    pr <- c["Precision"]
  }
  f1 <- 2 * ((re * pr) / (re + pr))
 
  
  if (flag_print){print(paste("MACRO F1 is: ", round(f1, 2)))}
  return(f1)
}



compute.A.conditional <- function(pred.matrix, i, j, ref.outcome) {
  # computes A(i|j), the probability that a randomly 
  # chosen member of class j has a lower estimated probability (or score) 
  # of belonging to class i than a randomly chosen member of class i
  
  # select predictions of class members
  i.idx <- which(ref.outcome == i)
  j.idx <- which(ref.outcome == j)
  pred.i <- pred.matrix[i.idx, i] # p(G = i) assigned to class i observations
  pred.j <- pred.matrix[j.idx, i] # p(G = i) assigned to class j observations
  all.preds <- c(pred.i, pred.j)
  classes <- c(rep(i, length(pred.i)), rep(j, length(pred.j)))
  o <- order(all.preds)
  classes.o <- classes[o]
  # Si: sum of ranks from class i observations
  Si <- sum(which(classes.o == i))
  ni <- length(i.idx)
  nj <- length(j.idx)
  # calculate A(i|j)
  A <- (Si - ((ni * (ni + 1))/2)) / (ni * nj)
  return(A)
}


# https://link.springer.com/article/10.1023/A:1010920819831
#Hand, D.J., Till, R.J. A Simple Generalisation of the Area Under the ROC Curve for Multiple Class Classification Problems. 
#Machine Learning 45, 171-186 (2001). https://doi.org/10.1023/A:1010920819831

multiclass.auc <- function(pred.matrix, ref.outcome) {
  labels <- colnames(pred.matrix)
  A.ij.cond <- utils::combn(labels, 2, function(x, pred.matrix, ref.outcome) {x
    i <- x[1]
    j <- x[2]
    A.ij <- compute.A.conditional(pred.matrix, i, j, ref.outcome)
    A.ji <- compute.A.conditional(pred.matrix, j, i, ref.outcome)
    pair <- paste0(i, "/", j)
    return(c(A.ij, A.ji))
  }, simplify = FALSE, pred.matrix = pred.matrix, ref.outcome = ref.outcome)
  c <- length(labels)
  pairs <- unlist(lapply(combn(labels, 2, simplify = FALSE), function(x) paste(x, collapse = "/")))
  A.mean <- unlist(lapply(A.ij.cond, mean))
  names(A.mean) <- pairs
  A.ij.joint <- sum(unlist(A.mean))
  M <- 2 / (c * (c-1)) * A.ij.joint 
  attr(M, "pair_AUCs") <- A.mean
  return(M)
}


# model <- NaiveBayes(iris.train$Species ~ ., data = iris.train[, -5])
# pred <- predict(model, iris.test[,-5], type='raw')
# pred.matrix <- pred$posterior
# ref.outcome <- iris.test$Species
# M <- multiclass.auc(pred.matrix, ref.outcome)
# print(paste0("Generalized AUC is: ", round(as.numeric(M), 3)))
# 
# print(attr(M, "pair_AUCs")) # pairwise AUCs