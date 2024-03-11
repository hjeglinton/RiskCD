
# remove.packages("reticulate")
#
# install.packages("reticulate", repos = "http://cran.us.r-project.org")
# library(reticulate)
#
# version <- "3.9.12"
# install_python(version)
# virtualenv_create(envname="FasterRisk-environment", version = version)
# use_virtualenv("FasterRisk-environment")
#
# py_install("fasterrisk", pip=TRUE, envname="FasterRisk-environment")

fasterrisk <- import("fasterrisk")
np <- import("numpy", convert=FALSE)

get_FR_deviance <- function(mod, y){

  # Get predicted probs and classes
  p <- mod$pred_test

  # Deviance
  p[p == 1] <- 0.99999
  p[p == 0] <- 0.00001
  dev <- -2*sum(y*log(p)+(1-y)*log(1-p))

  return(dev = dev)
}


#' Run FasterRisk
#'
#' Implements FasterRisk method to learn risk score model.
#' @param X_train Numeric matrix of training data coefficients. Must contain
#' an intercept column.
#' @param y_train Numeric vector of outcomes from the training data (0/1).
#' @param X_test Optional, numeric matrix of coefficients from the test data. If
#' not supplied, the function will return predicted probabilities using the
#' training data coefficients.
#' @return List with three attributes:
#' * `coef`: Numeric vector of logistic regression coefficients
#' * `integer_coef`: Numeric vector of integer coefficients.
#' * `pred_test`: Numeric vector with predicted probabilities calculated from the
#' `X_test` matrix, if supplied, or the `X_train` matrix.
run_FR <- function(X_train, y_train, X_test = NULL, lb, ub, k = dim(X_train)[2]-1) {

  if (is.null(X_test)) {
    X_test <- X_train
  }


  # Reformat outcome data to -1/1
  y_train <- case_when(y_train == 0 ~ -1,
                       y_train == 1 ~ 1)

  # Reformat to np arrays
  X_train <- np$array(X_train)
  y_train <- np$array(y_train)

  # Run FasterRisk method
  m <- fasterrisk$fasterrisk$RiskScoreOptimizer(X = X_train, y = y_train,
                                                k = k,
                                                lb = lb, ub = ub)
  m$optimize()


  # Get all solutions
  solutions = m$get_models()
  multipliers = solutions[1][[1]]
  sparseDiversePool_beta0_integer = solutions[2][[1]]
  sparseDiversePool_betas_integer = solutions[3][[1]]

  # Get top solution
  model_index = 1
  multiplier = multipliers[model_index]
  intercept = sparseDiversePool_beta0_integer[model_index]
  coefficients = np$array(sparseDiversePool_betas_integer[model_index, ])
  coefficients[0] = intercept

  # Calculate predicted probabilities of test set
  RiskScoreClassifier_m = fasterrisk$fasterrisk$RiskScoreClassifier(multiplier, intercept, coefficients)
  predicted_probs_train <- RiskScoreClassifier_m$predict_prob(X_train)
  predicted_probs_test <- RiskScoreClassifier_m$predict_prob(X_test)

  # Return list with logistic model, integer model, and predicted probabilities
  return(list(coef = as.numeric(coefficients/multiplier),
              integer_coef = as.numeric(coefficients),
              pred_train = predicted_probs_train,
              pred_test = predicted_probs_test))

}


run_FR_CV <- function(X, y, lb, ub, nfolds = 10,
                      num_k = min(floor(dim(X)[2]-1), 25), k_grid = NULL,
                      foldids = NULL, seed = NULL, parallel = FALSE) {
  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Get folds
  if (is.null(foldids) & is.null(nfolds)) stop("Must provide foldids or nfolds")
  if (is.null(foldids)){
    foldids <- sample(rep(seq(nfolds), length = nrow(X)))
  } else {
    nfolds <- max(foldids)
  }

  # Check at least 3 folds
  if (nfolds <= 3) stop("Must have more than 3 folds")

  # Check num_k less than number of columns
  if (num_k > dim(X)[2]-1) stop("num_k greater than number of columns")

  # Get k sequence
  #if (is.null(k_grid)){

  k_max <- dim(X)[2]-1
  k_min <- 0

  k_grid <- floor(seq(k_max, k_min, length.out=num_k))

  #}

  if (length(k_grid) < 1) stop("something wrong with k_grid")

  # Results data frame
  res_df <- data.frame(k = rep(k_grid, nfolds),
                       fold = sort(rep(1:nfolds, num_k)),
                       dev = rep(0, nfolds*num_k),
                       non_zeros = rep(0, nfolds*num_k))


  # Function to run for single fold and k
  fold_fcn <- function(k, foldid){

    X_train <- X[foldids != foldid, ]
    y_train <- y[foldids != foldid]
    X_test <- X[foldids == foldid,]
    y_test <- y[foldids == foldid]

    mod <- run_FR(X_train, y_train, X_test,
                  lb = lb, ub = ub, k = k)

    dev <- get_FR_deviance(mod, y_test)

    non_zeros <- sum(mod$integer_coef != 0)
    return(c(dev, non_zeros))
  }

  # Run through all folds
  # Parallel Programming. ! Must register parallel beforehand
  i = NULL # set global variable
  if (parallel) {

    outlist = foreach::foreach(i = 1:nrow(res_df)) %dopar%
      {
        fold_fcn(res_df[i,1],res_df[i,2])
      }
    res_df[,3:4] <- base::t(sapply(1:nrow(res_df), function(i) res_df[i,3:4] <- outlist[[i]]))
  } else {

    res_df[,3:4] <- base::t(sapply(1:nrow(res_df),
                                   function(i) fold_fcn(res_df$k[i],
                                                        res_df$fold[i])))
  }


  # Summarize
  dev = NULL # set global variables
  res_df_summary <- res_df %>%
    dplyr::group_by(k) %>%
    dplyr::summarize(mean_dev = mean(dev), sd_dev = stats::sd(dev),
                     mean_nonzeros = mean(non_zeros), sd_nonzeros = stats::sd(non_zeros))


  # Find k_min for deviance
  k_min_ind <- which.min(res_df_summary$mean_dev)
  k_min <- res_df_summary$k[k_min_ind]

  final <- list(results = res_df_summary, k_min = k_min)

  return(final)

}
