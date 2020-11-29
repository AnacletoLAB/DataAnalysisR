

mean_on_columns <- function(df_mat){
  return(apply(df_mat, 2, mean))
}


var_on_columns <- function(df_mat){
  var_on_col = apply(df_mat, 2, var)
#  se = sd/sqrt(nrow(df_mat))
  
  return(var_on_col)
}

se_on_columns <- function(df_mat){
  sd_on_col = apply(df_mat, 2, sd)
  se = sd_on_col/sqrt(nrow(df_mat))
  
  return(se)
}



min_max_norm <- function(x){
  #normalizzazione min max
  mins = apply(x, 2, min, na.rm = TRUE)
  maxs = apply(x, 2, max, na.rm = TRUE)
  for (nc in 1:ncol(x)){
    x[,nc] = (x[,nc]-mins[nc])/(maxs[nc]-mins[nc])
  }
  return(list("norm" = x, "min" = mins, "max" = maxs ))
}


norm_with_coeffs <- function(x, col_coeffs){
# normalizzo usando coefficienti dati
  min = col_coeffs[["min"]]
  max = col_coeffs[["max"]]
  for (nc in 1:ncol(x)){
    x[,nc] = (x[,nc]- col_coeffs[["min"]][nc]) /(col_coeffs[["max"]][nc]-col_coeffs[["min"]][nc])
  }
  return(x)
}