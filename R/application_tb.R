#devtools::install_github("hjeglinton/riskscores")

library(tidyverse)
library(riskscores)
setwd("/Users/hannaheglinton/Library/CloudStorage/OneDrive-BrownUniversity/Thesis/Experiments/R")
set.seed(1)

# Read in TB risk data
tb_risk_df <- read.csv("../data/tb/tb_risk_data.csv")

# Train/test split 


# Transform data into y and X matrices
y <- tb_risk_df[[1]]
X <- as.matrix(tb_risk_df[,-1])

# Assign stratified folds
foldids <- stratify_folds(y, nfolds = 10, seed = 1)

# Cross-validation to find lambda0
cv_results <- cv_risk_mod(X, y, foldids = foldids, nlambda = 100)
plot(cv_results)

# Fit model
mod <- risk_mod(X, y, lambda0 = cv_results$lambda_min)
mod$model_card
plot(mod, score_min = -50, score_max = 50)

get_metrics(mod)

# Get metrics on derivation data
y_pred <- predict(mod, type = "response") %>% as.vector()
roc_obj <- roc(y, y_pred)
threshold <- coords(roc_obj, "best",  ret = "threshold")
auc <- roc_obj$auc
brier <- mean((y_pred - y)^2)
precision <- coords(roc_obj, "best", ret = "precision")[[1]]
recall <- coords(roc_obj, "best", ret = "recall")[[1]]
f1_score <- 2 * (precision * recall) / (precision + recall)

brier <- mean((y_pred - y)^2)

# Get metrics on validation data (test/train split and/or use paper's derivation data)



# Get Baik predicted values 


