substitute_min_max <- function(df_in, min_vals, max_vals, string_to_print=''){
  df <- df_in
  col_names <- colnames(df)
  for (nf in 1:ncol(df)){
    name_col <- col_names[nf]
    cat(nf, ')', name_col, ', [', min_vals[nf], ', ',max_vals[nf], ']\n' )
    idx_min <-  (df[, nf] < min_vals[nf])
    if (any(idx_min)){
      cat(name_col, string_to_print, ' : values lower than observed min found and substituted.\n')
      cat( df[idx_min, nf], '\n')
      df[idx_min, nf] = min_vals[nf]
    }
    idx_max <- df[, nf] > max_vals[nf]
    if (any(idx_max)){
      cat(name_col, string_to_print, ' :values bigger than observed max found and substituted.\n')
      df[idx_max, nf] = max_vals[nf] 
    }
  }
  return(df)
}



extract_type <- function(data){
  
  library(stringr)
  
  cat('Extract data type from names', '\n')
  
  cnames = colnames(data)
  var_types <- data.frame(matrix(ncol = length(cnames), nrow = 1))
  colnames(var_types) = colnames(data)
  
  var_types[which(str_locate(str_to_lower(cnames), 'id')[,1]==1)] = 'ID'
  var_types[which(str_locate(str_to_lower(cnames), 'label')[,1]==1)] = cnames[which(str_locate(str_to_lower(cnames), 'label')[,1]>0)]
  var_types[which(str_locate(str_to_lower(cnames), 'ord')[,1]==1)] = 'ord'
  var_types[which(str_locate(str_to_lower(cnames), 'cat')[,1]==1)] = 'cat'
  var_types[which(str_locate(str_to_lower(cnames), 'int')[,1]==1)] = 'int'
  var_types[which(str_locate(str_to_lower(cnames), 'num')[,1]==1)] = 'num'
 
  return(var_types)
  
}


transform.data.with.type <- function(new_df){
  
  if (!is.data.frame(new_df)){
    new_df = as.data.frame(new_df)
  }
  cat('Convert data', '\n')
  
  data_types <- extract_type(new_df)
  
  nom_vars <- grep('cat', data_types[1,]) 
  ord_vars <- grep('ord', data_types[1,]) 
  int_vars <- grep('int', data_types[1,])
  double_vars <- grep('num', data_types[1, ])
  
  
  if (length(int_vars)>1){new_df[, int_vars] <- apply(apply(new_df[ , int_vars],2, as.numeric),2,round)}
  
  if (length(double_vars)>1){new_df[, double_vars] <- apply(apply(new_df[ , double_vars],2, as.numeric),2,round, 2)}
  
  if (length(nom_vars)>1){new_df[, nom_vars] <- apply(new_df[, nom_vars], 2, as.factor)}
  
  if (length(ord_vars)>1){ new_df[, ord_vars] <- apply(new_df[, ord_vars], 2, as.ordered)}
  
  return(new_df)
  
}


convert.rules <- function(ruleExec, new.idx){
   
  old_num = str_split(toString(1:length(new.idx)), pattern=', ', )[[1]]
  new_num = str_split(toString(new.idx), pattern=', ', )[[1]]
  cat('convert ', old_num, '  in:\n')
  cat(new_num, '\n')
  for (on in 1:length(old_num)){
    old_num[on] = paste(",",old_num[on],"]", sep="")
    new_num[on] = paste(",",new_num[on],"]", sep="")
    ruleExec[,1] = str_replace_all(ruleExec[,1], old_num[on], new_num[on])
  }
  return(ruleExec)
}


extract.rules.on.important <- function(ruleExec, keep.idx){
  keep_idx = str_split(toString(keep.idx), pattern=', ', )[[1]]
  cat('keep only:', '\n')
  cat(keep_idx)
  
  for (nk in 1:length(keep_idx)){
    if (nk ==1){keep = str_detect(ruleExec, paste(",",keep_idx[nk],"]", sep=""))}  
    else{ keep= keep | str_detect(ruleExec, paste(",",keep_idx[nk],"]", sep=""))}
  }
  return(which(keep))
}

convert_names <- function(learner.conditions, names.var){
  
  for (on in length(var.names):1){
   
    
    learner.conditions[,1] = 
        str_replace_all(learner.conditions[,1], 
                       paste(",",toString(on),"]",sep=""), var.names[on])
  }
 
  return(learner.conditions)
}


mean_on_cols <- function(df_mat){
  return(apply(df_mat, 2, mean))
}