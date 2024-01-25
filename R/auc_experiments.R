setwd("~/Documents/GitHub/thesis/R")

# Files in path
files <- list.files("../data/simulated/")

# Find data files
files_data <- files[grep("_data.csv", files)]
results <- data.frame(file = NA, auc = NA)
  
for (f in files_data) {
  # Read in data
  df <- read.csv(paste0("../data/simulated/",f))
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
  test_index <- createDataPartition(y, p = 0.3, list = FALSE)
  
  X_train <- X[-test_index,]
  X_test <- X[test_index,]
  y_train <- y[-test_index]
  y_test <- y[test_index]
  weights_train <- weights[-test_index]
  weights_test <- weights[test_index]
  
  
  # Logistic Regression
  df_glm <- data.frame(X_train, y = y_train)
  
  start_glm <- Sys.time()
  mod_glm <- glm(y ~ ., data = df_glm, weights = weights_train, family = "binomial")
  end_glm <- Sys.time()
  
  coef_glm <- mod_glm$coef %>% as.vector()
  pred_glm <- predict(mod_glm, data.frame(X_test), type = "response")
  
  roc_obj <- roc(y_test, pred_glm, quiet = TRUE)
  
  results <- bind_rows(results, data.frame(file = f, auc = roc_obj$auc[[1]])) 
  
}

results_auc <- separate(data = results, col = 1, 
                        into = c("sim", "n", "p", "prop_zero",
                                 "eps", "iter", "data"),
                        sep = "_") %>%
  group_by(n, p, prop_zero, eps) %>%
  summarize(auc = mean(auc)) %>%
  remove_missing()

View(results_auc %>%
  filter(n == 1000, p == 50))
