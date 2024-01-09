#install.packages("reticulate", repos = "http://cran.us.r-project.org")
#library(reticulate)

#version <- "3.9.12"
#install_python(version)
#virtualenv_create(envname="FasterRisk-environment", version = version)
#use_virtualenv("FasterRisk-environment")

#py_install("fasterrisk", pip=TRUE, envname="FasterRisk-environment")
fasterrisk <- import("fasterrisk")
np <- import("numpy", convert=FALSE)

FR_coef <- function(X, y, lb, ub) {
  
  y_train <- case_when(y == 0 ~ -1,
                       y == 1 ~ 1)
  
  X_train <- np$array(X)
  y_train <- np$array(y_train)
  
  m <- fasterrisk$fasterrisk$RiskScoreOptimizer(X = X_train, y = y_train, k = dim(X)[2]-1, 
                                                lb = lb, ub = ub)
  
  m$optimize()

  
  # get all top m solutions from the final diverse pool
  solutions = m$get_models()
  multipliers = solutions[1][[1]]
  sparseDiversePool_beta0_integer = solutions[2][[1]]
  sparseDiversePool_betas_integer = solutions[3][[1]]
  
  model_index = 1 # first model
  multiplier = multipliers[model_index]
  intercept = sparseDiversePool_beta0_integer[model_index]
  coefficients = np$array(sparseDiversePool_betas_integer[model_index, ])
  
  coefficients[0] = intercept

  return(as.numeric(coefficients/multiplier))
  
}


setwd("/Users/hannaheglinton/Library/CloudStorage/OneDrive-BrownUniversity/Thesis/Experiments/R")

# Read in TB risk data
tb_risk_df <- read.csv("../data/tb/tb_risk_data.csv")


# Transform data into y and X matrices
y <- tb_risk_df[[1]]
X <- as.matrix(tb_risk_df[,-1])

FR_coef(X, y, lb = -5, ub = 5)
