#devtools::install_github("hjeglinton/riskscores", build_vignettes = TRUE)

setwd("~/Documents/GitHub/thesis/R")
set.seed(1)

library(riskscores)
library(doParallel)
library(caret)
library(glmnet)
library(pROC)
library(tidyverse)

source("fasterrisk.R")



# Register parallel programming
registerDoParallel(cores=4)

#' Generate row for results dataframe
#' 
#' Summarizes model metrics into a dataframe row
#' @param betas Numeric vector of model coefficients, including the intercept. 
#' @param X Original dataframe (used for n and p values).
#' @param pred Numeric vector of predicted probabilities to compare to observed
#' outcomes. 
#' @param y Numeric vector of observed outcomes, of the same length as `pred`.
#' @param t1 Starting time.
#' @param t2 Ending time. 
#' @param file File name.
#' @param method Method name. 
#' @param lambda0 If applicable, lambda0 value used in model. 
#' @return Dataframe with a single row containing information on dataset, method,
#' and model metrics. 
get_result_row <- function(betas, X, pred, y, t1, t2, file, method, lambda0 = NA) {

  time_secs <- t2 - t1
  non_zeros <- sum(betas[-1] != 0, na.rm = TRUE)
  med_abs <- median(abs(betas[-1]), na.rm = TRUE)
  max_abs <- max(abs(betas[-1]), na.rm = TRUE)
  
  # Deviance
  pred[pred == 0] <- 1e-10
  pred[pred == 1] <- 1-1e-10
  dev <- -2*sum(y*log(pred)+(1-y)*log(1-pred))
  
  # Save ROC object 
  roc_obj <- roc(y, pred)
  
  precision <- coords(roc_obj, "best", ret = "precision")[[1]]
  recall <- coords(roc_obj, "best", ret = "recall")[[1]]
  
  # Return results
  return(data.frame(
    data = file,
    n = nrow(X), 
    p = ncol(X), 
    method = method,
    seconds = as.numeric(difftime(t2, t1, units = "secs")),
    lambda0 = lambda0,
    nonzeros = non_zeros,
    deviance = dev, 
    threshold = coords(roc_obj, "best",  ret = "threshold")[[1]],
    auc = roc_obj$auc[[1]],
    brier = mean((pred - as.numeric(y))^2),
    precision = precision,
    recall = recall,
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
                        method = character(), seconds = numeric(), 
                        lambda0 = numeric(), nonzeros = integer(), deviance = numeric(),
                        threshold = numeric(), auc = numeric(), brier = numeric(), 
                        precision = numeric(), recall = numeric(), 
                        f1_score = numeric(), accuracy = numeric(), 
                        sensitivity = numeric(), specificity = numeric(), 
                        med_abs = numeric(), max_abs = numeric())
  
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
    
    
    # NLLCD - no CV
    start_nllcd <- Sys.time()
    lambda0 <- 0
    mod_nllcd <- risk_mod(X_train, y_train, weights=weights_train, lambda0 = lambda0)
    coef_nllcd <- coef(mod_nllcd) %>% as.vector
    end_nllcd <- Sys.time()
    
    pred_nllcd <- predict(mod_nllcd, X_test, type = "response")[,1]
    res_nllcd <- get_result_row(coef_nllcd, X, pred_nllcd, y_test, start_nllcd,
                                end_nllcd, f, "NLLCD",
                                lambda0)
    
    # NLLCD - with CV
    start_nllcd_cv <- Sys.time()
    cv_results <- cv_risk_mod(X_train, y_train, weights=weights_train, 
      foldids = foldids, parallel = T)
    mod_nllcd_cv_min <- risk_mod(X_train, y_train, weights=weights_train, lambda0 = cv_results$lambda_min)
    coef_nllcd_cv_min <- coef(mod_nllcd_cv_min) %>% as.vector
    end_nllcd_cv <- Sys.time()
    
    mod_nllcd_cv_1se <- risk_mod(X_train, y_train, weights=weights_train, lambda0 = cv_results$lambda_1se)
    coef_nllcd_cv_1se <- coef(mod_nllcd_cv_1se) %>% as.vector
    
    pred_nllcd_cv_min <- predict(mod_nllcd_cv_min, X_test, type = "response")[,1]
    pred_nllcd_cv_1se <- predict(mod_nllcd_cv_1se, X_test, type = "response")[,1]
    
    
    res_nllcd_cv_min <- get_result_row(coef_nllcd_cv_min, X, pred_nllcd_cv_min, y_test, start_nllcd_cv,
                                end_nllcd_cv, f, "NLLCD with CV (lambda_min)",
                                cv_results$lambda_min)
    
    res_nllcd_cv_1se <- get_result_row(coef_nllcd_cv_1se, X, pred_nllcd_cv_1se, y_test, start_nllcd_cv,
                                       end_nllcd_cv, f, "NLLCD with CV (lambda_min)",
                                       cv_results$lambda_1se)
    
    # FasterRisk
    X_train_FR <- cbind(rep(1, nrow(X_train)), X_train)
    X_test_FR <- cbind(rep(1, nrow(X_test)), X_test)
    
    start_FR <- Sys.time()
    mod_FR <- run_FR(X_train_FR, y_train, X_test_FR, lb = -10, ub = 10) 
    end_FR <- Sys.time()
    
    res_FR <- get_result_row(mod_FR$integer_coef, X, mod_FR$pred_test, y_test, 
                             start_FR, end_FR, f, "FasterRisk", NA)
    
    # Lasso
    start_lasso <- Sys.time()
    mod_lasso <- cv.glmnet(X_train, y_train, foldids = foldids, 
                           weights = weights_train, alpha = 1, family = "binomial") 
    end_lasso <- Sys.time()
    
    coef_lasso <- as.vector(coef(mod_lasso, lambda = mod_lasso$lambda.min))
    pred_lasso <- predict(mod_lasso, X_test, "lambda.min", type = "response")[,1]
    
    res_lasso <- get_result_row(coef_lasso, X, pred_lasso, y_test,
                                start_lasso, end_lasso, f, "Lasso", 
                                mod_lasso$lambda.min)
    
    # Rounded Lasso 
    nonzero_lasso <- coef_lasso[coef_lasso != 0][-1]
    b0_lasso_scaled <- coef_lasso[1]/median(nonzero_lasso)
    coef_lasso_rounded <- c(b0_lasso_scaled, 
                            round(coef_lasso[-1]/median(nonzero_lasso), 0))
    scores <- as.vector(X_train_FR %*% coef_lasso_rounded)
    mod_scores <- glm(y_train ~ scores, family = "binomial")
    end_lasso_rounded <- Sys.time()
    
    scores_test <- data.frame(scores = X_test_FR %*% coef_lasso_rounded)
    pred_lasso_rounded <- predict(mod_scores, scores_test, type = "response")
    
    res_lasso_rounded <- get_result_row(coef_lasso_rounded, X, pred_lasso_rounded,
                                        y_test, start_lasso, end_lasso_rounded, 
                                        f, "Rounded Lasso", mod_lasso$lambda.min)

    
  
    # Logistic Regression
    df_glm <- data.frame(X_train, y = y_train)
    
    start_glm <- Sys.time()
    mod_glm <- glm(y ~ ., data = df_glm, weights = weights_train, family = "binomial")
    end_glm <- Sys.time()
    
    coef_glm <- mod_glm$coef %>% as.vector()
    pred_glm <- predict(mod_glm, data.frame(X_test), type = "response")
    
    res_glm <- get_result_row(coef_glm, X, pred_glm, y_test, 
                              start_glm, end_glm, f, "LR", NA)
    
    
    # Rounded Logistic Regression
    b0_glm_scaled <- coef_glm[1]/median(coef_glm[-1], na.rm = TRUE)
    coef_glm_rounded <- c(b0_glm_scaled, 
                          round(coef_glm[-1]/median(coef_glm[-1], na.rm = TRUE))) %>%
      replace(is.na(.), 0)
    
    scores <- as.vector(X_train_FR %*% coef_glm_rounded)
    mod_scores <- glm(y_train ~ scores, family = "binomial")
    end_glm_rounded <- Sys.time()
    
    scores_test <- data.frame(scores = X_test_FR %*% coef_glm_rounded)
    pred_glm_rounded <- predict(mod_scores, scores_test, type = "response")
    
    res_glm_rounded <- get_result_row(coef_glm_rounded, X, pred_glm_rounded,
                                        y_test, start_glm, end_glm_rounded, 
                                        f, "Rounded LR", NA)

    # Combine evaluation metrics
    results <- bind_rows(
      results, 
      res_FR,
      res_glm,
      res_glm_rounded,
      res_lasso,
      res_lasso_rounded,
      res_nllcd,
      res_nllcd_cv_min,
      res_nllcd_cv_1se
    )
  }
  write.csv(results, results_path, row.names=FALSE)
}

run_experiments(data_path = "../data/public/", results_path = "../results/public/results_public_withCV_2.csv")
