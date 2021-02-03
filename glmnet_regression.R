glmnet_regression <- function(data, LABEL, ext.folds, n_alpha = 51, nameLab = "LABEL"){
  # legge i dati nel file data_filename, i tipi nel file type_filename, i folds nella matrice folds
  # data_types contiene LABELS come tipo delle colonne da usare come labels
    library(glmnet)
    library(glmnetUtils)
    library(cutpointr)
    
    # provo tutti questi alpha
    if (n_alpha>1){
      alpha = seq(0, 1, len = n_alpha)^3
    }else{
      if (n_alpha ==0){ alpha = 0 }
      else{ if (n_alpha == 1){ alpha =1}else{alpha = n_alpha} }
    }
    
    
    cat('working on label', nameLab, '\n')
    
   
    num_classes = length(levels(LABEL))
    num_folds = length(ext.folds)   
  #     #vettore dei coefficienti/voti delle features computate da logistic regression with lasso contraint
  #     #Ogni riga contiene i voti ottenuti per quel fold sulle NF features
    coeffs <- matrix(0, nrow = num_folds, ncol = ncol(data)) 
    colnames(coeffs) <- colnames(data)
    rownames(coeffs) <-  as.character(1:num_folds)
     
     # tiene la conta dei voti
    # votes <- matrix(0, nrow = 1, ncol = ncol(data)) 
    # colnames(votes) <- colnames(data)
    # 
    alpha.min = matrix(data = NA, nrow = length(alpha), ncol = 1)
     
    
     
    results_glm = rep(levels(LABEL)[1], length = nrow(data))
    p.test = as.numeric(LABEL)-as.numeric(LABEL)
    # applico generalized linear regression (gaussian family) with mix of ridge-lasso constraint sui 10 folds
    for (num_f in 1:num_folds){
     
      cat('num external fold = ', num_f, '\n')
      ns = ext.folds[[num_f]]
      # # apro il file dei tipi in cui le features non sono ancora state selezionate da Boruta: così ho i nomi di tutte le features
      # # butto via le colonne TAG che non corrispondono a features
      data.train = data[-ns,]
      lab.train = LABEL[-ns]
      
      data.test = data[ns,]
      lab.test = LABEL[ns]
     
     
     # invece che applicare lasso o ridge provo n_alpha valori di alpha e prendo la coppia alpha e lambda che mi danno il minimo mse 
     
      all_vals_alpha <- 1:n_alpha
      all_cvfit_lasso <- cva.glmnet(y = lab.train, x = data.train,type.measure = "mse", nfolds = 10,  
                                    alpha = alpha, standardize = TRUE, family = "binomial")  
      
      # la struttura ritornata contiene, per ogni alpha, 
      #la auc ottenuta da tutti i valori testati per lambsa
      for (nalpha in 1:length(all_vals_alpha)){
        all_vals_alpha[nalpha] <- all_cvfit_lasso$modlist[[nalpha]]$cvm[which(all_cvfit_lasso$modlist[[nalpha]]$lambda == all_cvfit_lasso$modlist[[nalpha]]$lambda.min)]
      }
#      plot(alpha, all_vals_alpha)
      alpha.min[num_f] = alpha[which.min(all_vals_alpha)]
          # modello di regressione lineare con best_alpha
     cvfit_lasso = all_cvfit_lasso$modlist[[which.min(all_vals_alpha)]]
     
     all_vals_alpha = NULL
     
     obj.coeffs <- coef(cvfit_lasso, s = "lambda.min", alpha = alpha.min[num_f])
     
     # ora registro sia i coefficienti che i voti
     coeffs[num_f, obj.coeffs@i[2:length(obj.coeffs@i)]]<- obj.coeffs@x[2:length(obj.coeffs@x)]
     # aggiungo un punto togliendo l'intercetta
     # votes[obj.coeffs@i[2:length(obj.coeffs@i)]] <- votes[obj.coeffs@i[2:length(obj.coeffs@i)]]+1
     # 
     
     fit_train <- predict(cvfit_lasso, newx = data.train, s = "lambda.min")
     p.train = 1/(1+exp(-fit_train))
     
     fit_test <- predict(cvfit_lasso, newx = data.test, s = "lambda.min")
     p.test[ns] = 1/(1+exp(-fit_test))
     
     # ora devo predire i migliori cutoff tra classi
    
     opt_cut <- cutpointr(x = p.train, class = lab.train, boot_runs = 100, pos_class = levels(LABEL)[2],
                          boot_stratify = TRUE, method = maximize_gam_metric, metric = youden, criterion = "aicc")
       
    
     results_glm[ns[ which(p.test[ns]>=mean(opt_cut$boot[[1]]$optimal_cutpoint))]] = levels(LABEL)[2]
     
   
    }
    
 

    
       return(list("LABELS"= LABEL, "predictions" = results_glm, "probs" = p.test, "coeffs" = coeffs))
     
   }
  