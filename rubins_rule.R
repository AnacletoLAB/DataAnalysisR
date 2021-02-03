rubin_rule_pool <- function(imp){
#https://thomasleeper.com/Rcourse/Tutorials/mi.html

  # imp is a list with num_imputations elements
  # each element is a dataframe/matrix where rows are the repetitions of the classifier and each column is 
  # related to a different value (either performance value or importance of the feature)
  source(file.path('.', 'Utils_postprocess.R'))
  num_imputations = length(imp)
  
  # with the following function I get a matrix with 
  # num_imputation columns - one for each imputation run.
  # each row is related to one of the values in the imp dataframe/matrix 
  # (essentially, for each value I compute the mean over all the repetitions of the algorithm)
  mean_estimates <- sapply(imp, mean_on_columns)
  
  
  #Now I compute the mean over all the imputations
  grandm <- apply(mean_estimates,1,mean)
  grandm
  
#  now for each imputation, I compute the standard error of the classifiers' runs
  ses <- sapply(imp, var_on_columns)
  stderrs <- sapply(imp, se_on_columns)
 
  #To get the standard error of our multiple imputation estimate, 
  #we need to combine the standard errors of each of our estimates, 
  #so that estimates we need to start by getting the SEs of each imputed vector:
  #The within variance is the mean of the se for all the imputations   
  within <- apply(ses,1, mean)
  within_se <- apply(stderrs,1, mean)
  
    
  #To calculate the between-imputation VARIANCE, 
  #we calculate the sum of squared deviations of each imputed mean from the grand mean estimate:
  # FOR EACH IMPUTATION COMPUTINE THE SAMPLE VARIANCE
  between = within-within
  for (nv in 1:nrow(mean_estimates)) {
    between[nv] = sum((mean_estimates[nv,]-grandm[nv])^2)* (1/(num_imputations-1))
  }
  
  between_se = sqrt(between)/sqrt(num_imputations)
  #Then we sum the within- and between-imputation variances (multiply the latter by a small correction):
    
#  cat(between, '\n')
  
  grandvar <- within + (1+1/num_imputations)*between 
  grandse <- within_se + (1+1/num_imputations)*between_se
  return(list("mean" = grandm, "var" = grandvar, "se" = grandse, "between" = between, "within" = within))
}



mean_on_cols <- function(df_mat){
  return(apply(df_mat, 2, mean))
}