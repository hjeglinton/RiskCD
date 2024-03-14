#devtools::install_github("hjeglinton/riskscores", build_vignettes = TRUE)

setwd("~/Documents/GitHub/thesis/R")
set.seed(1)

library(riskscores)
library(doParallel)
library(caret)
library(glmnet)
library(pROC)
library(tidyverse)
library(parallel)

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
get_result_row_all_metrics <- function(betas, X, pred, y, t1, t2, file, method, lambda0 = NA) {

  time_secs <- t2 - t1
  non_zeros <- sum(betas[-1] != 0, na.rm = TRUE)
  med_abs <- median(abs(betas[-1]), na.rm = TRUE)
  max_abs <- max(abs(betas[-1]), na.rm = TRUE)

  # Deviance
  pred[pred == 0] <- 1e-10
  pred[pred == 1] <- 1-1e-10
  dev <- -2*sum(y*log(pred)+(1-y)*log(1-pred))

  # Save ROC object
  roc_train <- roc(y)
  roc_test <- roc(y, pred, quiet = TRUE)

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
#' @param results_path character path to folder in which to save results
#' @return results are saved as a csv file results_R.csv to results_path
run_experiments <- function(data_path, results_path, random_split = TRUE,
                            method = c("riskcd", "riskcd-cv",
                                       "fasterrisk", "fasterrisk-cv",
                                       "lasso", "lasso-rounded",
                                       "logistic", "logistic-rounded")){


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

    # NLLCD - no CV
    if ("riskcd" %in% method) {
      start_nllcd <- Sys.time()
      lambda0 <- 0
      mod_nllcd <- risk_mod(X_train, y_train, weights=weights_train, lambda0 = lambda0)
      coef_nllcd <- coef(mod_nllcd) %>% as.vector
      end_nllcd <- Sys.time()

      results <- bind_rows(results,
                           get_result_row(coef_nllcd, X_train, y_train,
                                  pred_train = predict(mod_nllcd, X_train, type = "response")[,1],
                                  X_test, y_test,
                                  pred_test = predict(mod_nllcd, X_test, type = "response")[,1],
                                  start_nllcd, end_nllcd, f, "NLLCD", lambda0))
    }


    # NLLCD - with CV
    if ("riskcd-cv" %in% method) {
      cv_results <- cv_risk_mod(X_train, y_train, weights=weights_train,
                                foldids = foldids, parallel = T)
      start_nllcd_cv <- Sys.time()
      mod_nllcd_cv_min <- risk_mod(X_train, y_train, weights=weights_train, lambda0 = cv_results$lambda_min)
      coef_nllcd_cv_min <- coef(mod_nllcd_cv_min) %>% as.vector
      end_nllcd_cv <- Sys.time()

      results <- bind_rows(results,
                           get_result_row(coef_nllcd_cv_min, X_train, y_train,
                                         pred_train = predict(mod_nllcd_cv_min, X_train, type = "response")[,1],
                                         X_test, y_test,
                                         pred_test = predict(mod_nllcd_cv_min, X_test, type = "response")[,1],
                                         start_nllcd_cv, end_nllcd_cv, f,
                                         "NLLCD with CV (lambda_min)", cv_results$lambda_min))
    }



    # FasterRisk - no CV
    if ("fasterrisk" %in% method) {


      start_FR <- Sys.time()
      mod_FR <- run_FR(X_train_FR, y_train, X_test_FR, lb = -10, ub = 10, k = dim(X_train_FR)[2]-1)
      end_FR <- Sys.time()

      results <- bind_rows(results,
                           get_result_row(mod_FR$integer_coef, X_train, y_train,
                               pred_train = mod_FR$pred_train,
                               X_test, y_test,
                               pred_test = mod_FR$pred_test,
                               start_FR, end_FR, f, "FasterRisk", NA))
    }

    # FasterRisk - with CV
    if ("fasterrisk-cv" %in% method) {


      FR_CV <- run_FR_CV(X_train_FR, y_train, -10, 10, foldids = foldids, parallel = FALSE, k_grid = NULL)
      start_FR_CV <- Sys.time()
      mod_FR_CV <- run_FR(X_train_FR, y_train, X_test_FR, lb = -10, ub = 10, k = FR_CV$k_min)
      end_FR_CV <- Sys.time()

      results <- bind_rows(results,
                           get_result_row(mod_FR_CV$integer_coef, X_train, y_train,
                                          pred_train = mod_FR_CV$pred_train,
                                          X_test, y_test,
                                          pred_test = mod_FR_CV$pred_test,
                                          start_FR_CV, end_FR_CV, f, "FasterRisk-CV", NA))
    }


    # Lasso
    if ("lasso" %in% method) {
      start_lasso <- Sys.time()
      mod_lasso <- cv.glmnet(X_train, y_train, foldids = foldids,
                             weights = weights_train, alpha = 1, family = "binomial")
      end_lasso <- Sys.time()

      coef_lasso <- as.vector(coef(mod_lasso, lambda = mod_lasso$lambda.min))

      results <- bind_rows(results,
                           get_result_row(coef_lasso, X_train, y_train,
                                  pred_train = predict(mod_lasso, X_train, "lambda.min", type = "response")[,1],
                                  X_test, y_test,
                                  pred_test = predict(mod_lasso, X_test, "lambda.min", type = "response")[,1],
                                  start_lasso, end_lasso, f, "Lasso",
                                  mod_lasso$lambda.min))
    }


    # Rounded Lasso
    if ("lasso-rounded" %in% method) {
      if (!("lasso" %in% method)) {
        start_lasso <- Sys.time()
        mod_lasso <- cv.glmnet(X_train, y_train, foldids = foldids,
                               weights = weights_train, alpha = 1, family = "binomial")
        coef_lasso <- as.vector(coef(mod_lasso, lambda = mod_lasso$lambda.min))
      }

      scalar <- ifelse(is_empty(coef_lasso[coef_lasso != 0][-1]),
                       1,
                       median(coef_lasso[coef_lasso != 0][-1], na.rm = TRUE))
      b0_lasso_scaled <- coef_lasso[1]/scalar
      coef_lasso_rounded <- c(b0_lasso_scaled,
                              round(coef_lasso[-1]/scalar, 0))
      scores <- as.vector(X_train_FR %*% coef_lasso_rounded)
      mod_scores <- glm(y_train ~ scores, family = "binomial")
      end_lasso_rounded <- Sys.time()

      scores_test <- data.frame(scores = X_test_FR %*% coef_lasso_rounded)

      results <- bind_rows(results,
                           get_result_row(coef_lasso_rounded, X_train, y_train,
                                          pred_train = predict(mod_scores, type = "response"),
                                          X_test, y_test,
                                          pred_test = predict(mod_scores, scores_test, type = "response"),
                                          start_lasso, end_lasso_rounded,
                                          f, "Rounded Lasso", mod_lasso$lambda.min))
    }


    # Logistic Regression
    if ("logistic" %in% method) {
      df_glm <- data.frame(X_train, y = y_train)

      start_glm <- Sys.time()
      mod_glm <- glm(y ~ ., data = df_glm, weights = weights_train, family = "binomial")
      end_glm <- Sys.time()

      coef_glm <- mod_glm$coef %>% as.vector()


      results <- bind_rows(results,
                           get_result_row(coef_glm, X_train, y_train,
                                pred_train = predict(mod_glm, type = "response"),
                                X_test, y_test,
                                pred_test = predict(mod_glm, data.frame(X_test), type = "response"),
                                start_glm, end_glm, f, "LR", NA))
    }



    # Rounded Logistic Regression
    if ("logistic-rounded" %in% method) {
      if (!("logistic") %in% method) {
        start_glm <- Sys.time()
        mod_glm <- glm(y ~ ., data = df_glm, weights = weights_train, family = "binomial")
        coef_glm <- mod_glm$coef %>% as.vector()
      }
      b0_glm_scaled <- coef_glm[1]/median(coef_glm[-1], na.rm = TRUE)
      coef_glm_rounded <- c(b0_glm_scaled,
                            round(coef_glm[-1]/median(coef_glm[-1], na.rm = TRUE))) %>%
        replace(is.na(.), 0)

      scores <- as.vector(X_train_FR %*% coef_glm_rounded)
      mod_scores <- glm(y_train ~ scores, family = "binomial")
      end_glm_rounded <- Sys.time()

      scores_test <- data.frame(scores = X_test_FR %*% coef_glm_rounded)

      results <- bind_rows(results,
                           get_result_row(coef_glm_rounded, X_train, y_train,
                                          pred_train = predict(mod_scores, type = "response"),
                                          X_test, y_test,
                                          pred_test = predict(mod_scores, scores_test, type = "response"),
                                          start_glm, end_glm_rounded,
                                          f, "Rounded LR", NA))

    }

    return(results %>% filter_all(any_vars(!is.na(.))))
  }, mc.cores = 4)


  combine_results <- do.call(rbind, lapply(res, data.frame))

  write.csv(combine_results, results_path, row.names=FALSE)
}


# Public datasets
#run_experiments(data_path = "../data/public/", results_path = "../results/public/results_public.csv", random_split = TRUE)


# Simulated data

run_experiments(data_path = "../data/simulated/",
                results_path = "../results/simulated/results_sim_0314_FR_noCV.csv",
                random_split = FALSE,
                method = c("fasterrisk"))









