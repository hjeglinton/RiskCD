setwd("~/Documents/GitHub/thesis/R")

library(riskscores)
library(doParallel)
library(caret)
library(glmnet)


# Register parallel programming
registerDoParallel(cores=4)


get_result_row <- function(betas, X, pred, y, t1, t2, file, method, lambda0) {

  time_secs <- t2 - t1
  non_zeros <- sum(betas[-1] != 0)
  med_abs <- median(abs(mod$beta[-1]))
  max_abs <- max(abs(mod$beta[-1]))
  
  # Deviance
  dev <- -2*sum(y*log(pred)+(1-y)*log(1-pred))
  
  # Save ROC object 
  roc_obj <- roc(y, pred)
  
  # Return results
  return(data.frame(
    data = file,
    n = nrow(X), 
    p = ncol(X) - 1, 
    method = method,
    seconds = as.numeric(difftime(t2, t1, units = "secs")),
    lambda0 = lambda0,
    nonzeros = non_zeros,
    threshold = coords(roc_obj, "best",  ret = "threshold")[[1]],
    auc = roc_obj$auc[[1]],
    brier = mean((pred - as.numeric(y) + 1)^2),
    precision = coords(roc_obj, "best", ret = "precision")[[1]],
    recall = coords(roc_obj, "best", ret = "recall")[[1]],
    f1_score = 2 * (precision * recall) / (precision + recall),
    accuracy = coords(roc_obj, "best", ret = "accuracy")[[1]],
    sensitivity = coords(roc_obj, "best", ret = "sensitivity")[[1]],
    specificity = coords(roc_obj, "best", ret = "specificity")[[1]],
    med_abs = med_abs,
    max_abs = max_abs
  ))
  
} 

#' Risk model algorithm experiments
#' 
#' Iterates through all csv files in the path and runs risk mod
#' @param data_path character path to folder of csv files (assumes
#' files are stored as (x)_data.csv and (x)_weights.csv)
#' @param results_path character path to folder in which to save results
#' @return results are saved as a csv file results_R.csv to results_path
run_experiments <- function(data_path, results_path){

  # Files in path
  files <- list.files(data_path)
  results <- data.frame(data = character(), n = numeric(), p = numeric(),
                        lambda0 = numeric(), non_zeros = numeric(), 
                        med_abs = numeric(), max_abs = numeric(),
                        acc = numeric(), sens = numeric(), spec = numeric(),
                        sec = numeric())
  
  # Iterate through files
  for (f in files){
    if (length(grep("_data.csv", f)) == 0) next
    
    # Print for ease
    print(paste0(data_path,f))
    
    # Read in data
    df <- read.csv(paste0(data_path,f))
    y <- df[[1]]
    X <- as.matrix(df[,2:ncol(df)])
    #X <- cbind(rep(1,nrow(X)), X) # adds intercept column
    
    # Add weights file if needed
    weights <- rep(1, nrow(X))
    weights_file <- paste0(substr(f,1,nchar(f)-8),"_weights.csv")
    if (file.exists(weights_file)){
      weights <- read.csv(weights_file)
      weights <- weights[[1]]
    }
    
    # Test train split
    test_index <- createDataPartition(y, p = 0.3, list = FALSE)
    
    X_train <- X[-test_index,]
    X_test <- X[test_index,]
    y_train <- y[-test_index]
    y_test <- y[test_index]
    weights_train <- weights[-test_index]
    weights_test <- weights[test_index]
    
    # Stratify folds 
    foldids <- stratify_folds(y_train, nfolds = 5, seed = 1)
    
    
    # NLLCD
    start_nllcd <- Sys.time()
    lambda0 <- cv_risk_mod(X_train, y_train, weights=weights_train, 
                           foldids = foldids, parallel = T)$lambda_min
    #lambda0 <- 0
    mod_nllcd <- risk_mod(X_train, y_train, weights=weights_train, lambda0 = lambda0)
    coef_nllcd <- coef(mod_nllcd$glm_mod) %>% as.vector
    stop_nllcd <- Sys.time()
    
    pred_nllcd <- predict(mod_nllcd, X_test, type = "response")[,1]
    res_nllcd <- get_result_row(coef_nllcd, X, pred_nllcd, y_test, s1, e1, f, "NLLCD",
                                lambda0)
    
    # FasterRisk
    source("fasterrisk.R")
    s2 <- Sys.time()
    X_train_FR <- cbind(rep(1, nrow(X_train)), X_train)
    coef_FR <- FR_coef(X_train_FR, y_train, lb = -10, ub = 10) %>% as.vector()
    e2 <- Sys.time()
    
    
    
    # Lasso
    s3 <- Sys.time()
    lasso_mod <- cv.glmnet(X_train, y_train, foldids = foldids, 
                           weights = weights_train, alpha = 1, family = "binomial") 
    coef_lasso <- as.vector(coef(lasso_mod, lambda = lasso_mod$lambda.min))
    e3 <- Sys.time()
    
    # Rounded Lasso 
    nonzero_lasso <- coef_lasso[coef_lasso != 0][-1]
    b0_lasso_scaled <- coef_lasso[1]/median(nonzero_lasso)
    coef_lasso_rounded <- c(b0_lasso_scaled, 
                            round(coef_lasso[-1]/median(nonzero_lasso), 0))
    scores <- X_train_FR %*% coef_lasso_rounded
    mod_scores <- glm(y_train ~ scores, family = "binomial")$coefficients


    # Logistic Regression
    s4 <- Sys.time()
    coef_glm <- glm(y_train ~ X_train, weights = weights_train, family = "binomial")$coef %>%
      as.vector()
    
    b0_glm_scaled <- coef_glm[1]/median(coef_glm)
    coef_glm_rounded <- c(b0_glm_scaled, 
                          round(coef_glm[-1]/median(coef_glm)))
    e4 <- Sys.time()
    

    
    # Get evaluation metrics
    results <- bind_rows(
      get_result_row()
    )
    
    
    
    
  }
  write.csv(results, paste0(results_path, "results_cv_R.csv"), row.names=FALSE)
}

run_experiments(data_path = "../data/public/", results_path = "../results/public/")
