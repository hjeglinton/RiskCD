#devtools::install_github("hjeglinton/riskscores", build_vignettes = TRUE)

setwd("/Users/Alice/Dropbox/RiskCD/R")
# setwd("/oscar/home/stang52/riskscores/R")

source("fasterrisk.R")
set.seed(1)

library(riskscores)
library(doParallel)
library(caret)
library(glmnet)
library(pROC)
library(tidyverse)
library(parallel)

# Register parallel programming
registerDoParallel(cores=4)



#' Generate row for results dataframe
#'
#' Summarizes model metrics into a dataframe row
#' @param betas Numeric vector of model coefficients, including the intercept.
#' @param X_train Matrix of training covariates.
#' @param y_train Vector of training outcomes.
#' @param pred_train Numeric vector of predicted probabilities (training set).
#' @param X_test Matrix of test covariates.
#' @param y_test Vector of test outcomes.
#' @param pred_test Numeric vector of predicted probabilities (test set).
#' @param t1 Starting time.
#' @param t2 Ending time.
#' @param folder File name.
#' @param method Method name.
#' @param lambda0 If applicable, lambda0 value used in model.
#' @return Dataframe with a single row containing information on dataset, method,
#' and model metrics.
get_result_row <- function(betas, X_train, y_train, pred_train,
                           X_test, y_test, pred_test, t1, t2,
                           file, method, lambda0 = NA) {
  
  time_secs <- t2 - t1
  num_nonzeros <- sum(betas[-1] != 0, na.rm = TRUE)
  which_nonzeros <- which(betas[-1] != 0)
  med_abs <- median(abs(betas[-1]), na.rm = TRUE)
  max_abs <- max(abs(betas[-1]), na.rm = TRUE)
  
  # Save ROC objects
  roc_train <- roc(y_train, pred_train, quiet = TRUE)
  roc_test <- roc(y_test, pred_test, quiet = TRUE)
  
  # Return results
  return(data.frame(
    data = file,
    n = nrow(X_train),
    p = ncol(X_train),
    method = method,
    seconds = as.numeric(difftime(t2, t1, units = "secs")),
    lambda0 = lambda0,
    num_nonzeros = num_nonzeros,
    which_nonzeros = paste(which_nonzeros, collapse = "_"),
    auc_train = roc_train$auc[[1]],
    auc_test = roc_test$auc[[1]],
    brier_train = mean((pred_train - as.numeric(y_train))^2),
    brier_test = mean((pred_test - as.numeric(y_test))^2),
    med_abs = med_abs,
    max_abs = max_abs
  ))
  
}

#' Risk model algorithm experiments
#'
#' Iterates through all csv files in the path and runs risk mod
#' @param data_path character path to folder of csv files (assumes
#' files are stored as (x)_data.csv and (x)_weights.csv)
#' @param folder character path to folder in which to save results
#' @return results are saved as a csv file
run_experiments <- function(data_path, results_path, random_split = TRUE,
                            method = c("riskcd", "riskcd-cv",
                                       "annealscore", "annealscore-cv",
                                       "fasterrisk", "fasterrisk-cv",
                                       "lasso", "lasso-rounded",
                                       "logistic", "logistic-rounded"),
                            a = -10, b = 10){
  

  # Files in path
  files <- list.files(data_path)
  
  # Find data files
  files_data <- files[grep("_data.csv", files)]
  
  res <- mclapply(files_data, function(f) {

    # Initialize results
    results <- data.frame(data = NA,
                          n = NA,
                          p = NA,
                          method = NA,
                          seconds = NA,
                          lambda0 = NA,
                          num_nonzeros = NA,
                          which_nonzeros = NA,
                          auc_train = NA,
                          auc_test = NA,
                          brier_train = NA,
                          brier_test = NA,
                          med_abs = NA,
                          max_abs = NA)
    
    # Read in data
    df <- read.csv(paste0(data_path,f))
    
    if (colnames(df)[1] == "train") {
      train <- df$train
      df <- df %>% select(-train)
    }
    
    y <- df[[1]]
    X <- as.matrix(df[,2:ncol(df)])
    
    # Add weights file if needed
    weights <- rep(1, nrow(X))
    weights_file <- paste0(substr(f,1,nchar(f)-8),"_weights.csv")
    if (file.exists(weights_file)){
      weights <- read.csv(weights_file)
      weights <- weights[[1]]
    }
    
    # Test train split
    if (random_split == TRUE) {
      set.seed(2)
      test_index <- createDataPartition(y, p = 0.3, list = FALSE)
    } else {
      test_index <- which(train == 0)
    }
    
    X_train <- X[-test_index,]
    X_test <- X[test_index,]
    y_train <- y[-test_index]
    y_test <- y[test_index]
    weights_train <- weights[-test_index]
    weights_test <- weights[test_index]
    
    X_train_FR <- cbind(rep(1, nrow(X_train)), X_train)
    X_test_FR <- cbind(rep(1, nrow(X_test)), X_test)
    
    # Stratify folds
    foldids <- stratify_folds(y_train, nfolds = 5, seed = 1)
    
    # Simulated annealing
    if ("annealscore" %in% method) {
      start_annealscore <- Sys.time()
      mod_annealscore <- risk_mod(X_train, y_train, weights=weights_train, 
                                  method="annealscore", a = a, b=b)
      coef_annealscore <- coef(mod_annealscore) %>% as.vector
      end_annealscore <- Sys.time()
      
      row_annealscore <- get_result_row(coef_annealscore, X_train, y_train,
                                  pred_train = predict(mod_annealscore, X_train, type = "response")[,1],
                                  X_test, y_test,
                                  pred_test = predict(mod_annealscore, X_test, type = "response")[,1],
                                  start_annealscore, end_annealscore, f, "AnnealScore", 0)
      
      results <- tryCatch(bind_rows(results, row_annealscore), error = results)
    }
    
    # Simulated annealing -cv
    if ("annealscore-cv" %in% method) {
      start_annealscore_cv <- Sys.time()
      cv_results <- cv_risk_mod(X_train, y_train, weights=weights_train,
                                   foldids = foldids, parallel = T, method="annealscore",
                                a = a, b = b)
      mod_annealscore_cv_min <- risk_mod(X_train, y_train, weights=weights_train, 
                                         lambda0 = cv_results$lambda_min, 
                                         method="annealscore", a = a, b =b)
      coef_annealscore_cv_min <- coef(mod_annealscore_cv_min) %>% as.vector
      end_annealscore_cv <- Sys.time()

      row_annealscore_cv <- get_result_row(coef_annealscore_cv_min, X_train, y_train,
                                      pred_train = predict(mod_annealscore_cv_min, X_train, type = "response")[,1],
                                      X_test, y_test,
                                      pred_test = predict(mod_annealscore_cv_min, X_test, type = "response")[,1],
                                      start_annealscore_cv, end_annealscore_cv, f,
                                      "AnnealScore-CV", cv_results$lambda_min)

      results <- tryCatch(bind_rows(results, row_annealscore_cv), error = results)
    }
    
    # RiskCD - no CV
    if ("riskcd" %in% method) {
      start_nllcd <- Sys.time()
      mod_nllcd <- risk_mod(X_train, y_train, weights=weights_train, 
                            method="riskcd", a = a, b =b)
      coef_nllcd <- coef(mod_nllcd) %>% as.vector
      end_nllcd <- Sys.time()
      
      row_nllcd <- get_result_row(coef_nllcd, X_train, y_train,
                                  pred_train = predict(mod_nllcd, X_train, type = "response")[,1],
                                  X_test, y_test,
                                  pred_test = predict(mod_nllcd, X_test, type = "response")[,1],
                                  start_nllcd, end_nllcd, f, "RiskCD", 0)
      
      results <- tryCatch(bind_rows(results, row_nllcd), error = results)
    }
    
    
    # RiskCD - with CV
    if ("riskcd-cv" %in% method) {
      start_nllcd_cv <- Sys.time()
      cv_results <- cv_risk_mod(X_train, y_train, weights=weights_train,
                                foldids = foldids, parallel = T, method="riskcd",
                                a = a, b = b)
      mod_nllcd_cv_min <- risk_mod(X_train, y_train, weights=weights_train, 
                                   lambda0 = cv_results$lambda_min, 
                                   method="riskcd", a = a, b = b)
      coef_nllcd_cv_min <- coef(mod_nllcd_cv_min) %>% as.vector
      end_nllcd_cv <- Sys.time()
      
      row_nllcd_cv <- get_result_row(coef_nllcd_cv_min, X_train, y_train,
                                     pred_train = predict(mod_nllcd_cv_min, X_train, type = "response")[,1],
                                     X_test, y_test,
                                     pred_test = predict(mod_nllcd_cv_min, X_test, type = "response")[,1],
                                     start_nllcd_cv, end_nllcd_cv, f,
                                     "RiskCD-CV", cv_results$lambda_min)
      
      results <- tryCatch(bind_rows(results, row_nllcd_cv), error = results)
      
    }
    
    
    
    
    # FasterRisk - no CV
    if ("fasterrisk" %in% method) {
      
      start_FR <- Sys.time()
      mod_FR <- run_FR(X_train_FR, y_train, X_test_FR, lb = a, ub = b, k = dim(X_train_FR)[2]-1)
      end_FR <- Sys.time()
      
      row_FR <- get_result_row(mod_FR$integer_coef, X_train, y_train,
                               pred_train = mod_FR$pred_train,
                               X_test, y_test,
                               pred_test = mod_FR$pred_test,
                               start_FR, end_FR, f, "FasterRisk", NA)
      
      results <- tryCatch(bind_rows(results, row_FR), error = results)
    }
    
    # FasterRisk - with CV
    if ("fasterrisk-cv" %in% method) {
      
      if (ncol(X_train_FR) < 50) {
        
        start_FR_CV <- Sys.time()
        FR_CV <- run_FR_CV(X_train_FR, y_train, a, b, foldids = foldids,
                           parallel = FALSE, k_grid = seq(5, ncol(X_train_FR)-1, 5))
        mod_FR_CV <- run_FR(X_train_FR, y_train, X_test_FR, lb = a, ub = b, k = FR_CV$k_min)
        end_FR_CV <- Sys.time()
        
        row_FR_CV <- get_result_row(mod_FR_CV$integer_coef, X_train, y_train,
                                    pred_train = mod_FR_CV$pred_train,
                                    X_test, y_test,
                                    pred_test = mod_FR_CV$pred_test,
                                    start_FR_CV, end_FR_CV, f, "FasterRisk-CV", FR_CV$k_min)
        
        results <- tryCatch(bind_rows(results, row_FR_CV), error = results)
        
      }
      
      
      
    }
    
    
    
    # Lasso
    if ("lasso" %in% method) {
      start_lasso <- Sys.time()
      mod_lasso <- cv.glmnet(X_train, y_train, foldids = foldids,
                             weights = weights_train, alpha = 1, family = "binomial")
      end_lasso <- Sys.time()
      
      coef_lasso <- as.vector(coef(mod_lasso, lambda = mod_lasso$lambda.min))
      
      row_lasso <- get_result_row(coef_lasso, X_train, y_train,
                                  pred_train = predict(mod_lasso, X_train, "lambda.min", type = "response")[,1],
                                  X_test, y_test,
                                  pred_test = predict(mod_lasso, X_test, "lambda.min", type = "response")[,1],
                                  start_lasso, end_lasso, f, "Lasso",
                                  mod_lasso$lambda.min)
      
      results <- tryCatch(bind_rows(results, row_lasso), error = results)
    }
    
    
    # Rounded Lasso
    if ("lasso-rounded" %in% method) {
      if (!("lasso" %in% method)) {
        start_lasso <- Sys.time()
        mod_lasso <- cv.glmnet(X_train, y_train, foldids = foldids,
                               weights = weights_train, alpha = 1, family = "binomial")
        coef_lasso <- as.vector(coef(mod_lasso, lambda = mod_lasso$lambda.min))
      }
      
      nonzero_lasso <- coef_lasso[coef_lasso != 0][-1]
      
      scalar <- ifelse(is_empty(nonzero_lasso),
                       1,
                       max(abs(nonzero_lasso))/max(c(-a, b)))
      b0_lasso_scaled <- coef_lasso[1]/scalar
      coef_lasso_rounded <- c(b0_lasso_scaled,
                              round(coef_lasso[-1]/scalar, 0))
      scores <- as.vector(X_train_FR %*% coef_lasso_rounded)
      mod_scores <- glm(y_train ~ scores, family = "binomial")
      end_lasso_rounded <- Sys.time()
      
      scores_test <- data.frame(scores = X_test_FR %*% coef_lasso_rounded)
      
      
      row_lasso_rounded <- get_result_row(coef_lasso_rounded, X_train, y_train,
                                          pred_train = predict(mod_scores, type = "response"),
                                          X_test, y_test,
                                          pred_test = predict(mod_scores, scores_test, type = "response"),
                                          start_lasso, end_lasso_rounded,
                                          f, "Rounded Lasso", mod_lasso$lambda.min)
      
      results <- tryCatch(bind_rows(results, row_lasso_rounded), error = results)
    }
    
    
    # Logistic Regression
    if ("logistic" %in% method) {
      df_glm <- data.frame(X_train, y = y_train)
      
      start_glm <- Sys.time()
      mod_glm <- glm(y ~ ., data = df_glm, weights = weights_train, family = "binomial")
      end_glm <- Sys.time()
      
      coef_glm <- mod_glm$coef %>% as.vector()
      
      
      row_glm <- get_result_row(coef_glm, X_train, y_train,
                                pred_train = predict(mod_glm, type = "response"),
                                X_test, y_test,
                                pred_test = predict(mod_glm, data.frame(X_test), type = "response"),
                                start_glm, end_glm, f, "LR", NA)
      
      results <- tryCatch(bind_rows(results, row_glm), error = results)
    }
    
    
    
    # Rounded Logistic Regression
    if ("logistic-rounded" %in% method) {
      if (!("logistic") %in% method) {
        df_glm <- data.frame(X_train, y = y_train)
        start_glm <- Sys.time()
        mod_glm <- glm(y ~ ., data = df_glm, weights = weights_train, family = "binomial")
        coef_glm <- mod_glm$coef %>% as.vector()
      }
      
      scalar <- max(abs(coef_glm[-1]), na.rm = TRUE)/max(c(-a, b))
      coef_glm_rounded <- c(coef_glm[1]/scalar,
                            round(coef_glm[-1]/scalar, 0)) %>%
        replace(is.na(.), 0)
      
      scores <- as.vector(X_train_FR %*% coef_glm_rounded)
      mod_scores <- glm(y_train ~ scores, family = "binomial")
      end_glm_rounded <- Sys.time()
      
      scores_test <- data.frame(scores = X_test_FR %*% coef_glm_rounded)
      
      row_glm_rounded <- get_result_row(coef_glm_rounded, X_train, y_train,
                                        pred_train = predict(mod_scores, type = "response"),
                                        X_test, y_test,
                                        pred_test = predict(mod_scores, scores_test, type = "response"),
                                        start_glm, end_glm_rounded,
                                        f, "Rounded LR", NA)
      
      results <- tryCatch(bind_rows(results, row_glm_rounded), error = results)
      
      
    }
    
    return(results %>% filter_all(any_vars(!is.na(.))))
  }, mc.cores = 4)
  
  
   combine_results <- do.call(rbind, lapply(res, data.frame))

  write.csv(combine_results, results_path, row.names=FALSE)
}

# Simulated data
run_experiments(data_path = "../data/simulated-new/",
                results_path = "../results-new/results_noCV_10.csv",
                random_split = FALSE,
                method = c("riskcd",
                           "annealscore",
                           "lasso", "lasso-rounded",
                           "logistic", "logistic-rounded",
                           "fasterrisk"),
                a = -10, b = 10)
run_experiments(data_path = "../data/simulated-new/",
                results_path = "../results-new/results_noCV_5.csv",
                random_split = FALSE,
                method = c("riskcd",
                           "annealscore",
                           "lasso", "lasso-rounded",
                           "logistic", "logistic-rounded",
                           "fasterrisk"),
                a = -5, b = 5)
run_experiments(data_path = "../data/simulated-new/",
                results_path = "../results-new/results_noCV_3.csv",
                random_split = FALSE,
                method = c("riskcd",
                           "annealscore",
                           "lasso", "lasso-rounded",
                           "logistic", "logistic-rounded",
                           "fasterrisk"),
                a = 0, b = 3)
# run_experiments(data_path = "../data/",
#                 results_path = "../results/manuscript/results_fasterriskCV.csv",
#                 random_split = FALSE,
#                 method = c("fasterrisk-cv"))
run_experiments(data_path = "../data/simulated-new/",
                 results_path = "../results-new/results_CV.csv",
                 random_split = FALSE,
                 method = c("riskcd-cv","annealscore-cv"),
                a = -5, b=5)
run_experiments(data_path = "../data/simulated-new/",
                                 results_path = "../results-new/results_fasterriskCV.csv",
                                 random_split = FALSE,
                                 method = c("fasterrisk-cv"),
                a = -5, b=5)
# run_experiments(data_path = "../data/",
#                results_path = "../results/manuscript/results_fasterriskCV_p50.csv",
#                random_split = FALSE,
#                method = c("fasterrisk-cv"))
