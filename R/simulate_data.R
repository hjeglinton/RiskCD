library(tidyverse)
#library(MASS) # for multinomal 


#' Simulate data and coefficients
#' 
#' @param n Number of observations.
#' @param p Total number of variables.
#' @param p_zero Proportion of variables with a coefficient of zero.
#' @param eps Value between 0 and 1 representing amount of noise.
#' @return List with simulated covariates (x), outcome (y), coefficients (coef),
#' and intercept. 
sim_data <- function(n, p, prop_zero, snr){
  
  # Simulate coefficients
  p_zero <- floor(p*prop_zero)
  p_nonzero <- p - p_zero
  
  coef_small <- runif(floor(p_nonzero*(3/4)), 0.2, 1.5)
  coef_large <- runif(ceiling(p_nonzero*(1/4)), 2, 5)
  coef <- c(coef_small, coef_large, rep(0, p_zero))
  
  # Simulate covariates
  x <- matrix(0, nrow = n + (10*p), ncol = p)
  for (i in 1:p){
    x[,i] <- rbinom(n + (10*p), 1, runif(1, 0.1, 0.9))
  }
  
  # Get intercept
  vals <- x %*% coef
  intercept <- -mean(vals)
  
  # Get epsilon values
  sd_epsilon <- sqrt(var(vals, na.rm = TRUE) / snr)
  eps <- matrix(rnorm(n + (10*p), mean = 0, sd = sd_epsilon), ncol = 1)
  
  # Simulate outcome
  vals <- vals + intercept + eps
  probs <- exp(vals)/(1 + exp(vals))
  y <- rbinom(n + 10*p, 1, prob = probs)
  
  return(list(x = x, y = y,coef = coef, intercept = intercept))
}


#' Generate data and saves two files (one with data, one with coefficients)
#' 
#' @param n Number of observations.
#' @param p Total number of variables.
#' @param prop_zero Proportion of variables with a coefficient of zero.
#' @param eps Value between 0 and 1 representing amount of noise.
#' @param folder Path to folder in which to save files (include ending slash). 
#' @param filename File name prefix.
#' @return Saves two csv files with _data.csv and _coef.csv suffixes. 
gen_data <- function(n, p, prop_zero, snr, folder, filename){

  # Simulate data
  data <- sim_data(n, p, prop_zero, snr)
  
  # data file
  df <- as.data.frame(data$x)
  df$y <- data$y
  df$train <- c(rep(1, n), rep(0, 10*p))
  df <- df %>%
    select(train, y, everything())
  write.csv(df, paste0(folder, filename,"_data.csv"), row.names=FALSE)
  
  # coefficients file
  coef_df <- data.frame(names = c("Intercept", names(df)[3:length(df)]),
                        vals = c(data$intercept, data$coef))
  write.csv(coef_df, paste0(folder,filename,"_coef.csv"), row.names=FALSE)
  
}

# Example of generating data and saving
set.seed(5)
n = c(100, 500, 1000, 5000)
p = c(10, 25, 50)
prop_zero = c(0, 0.5)
snr = c(1, 3, 6)
n_iter = 10

for (j in 1:length(n)){
  for (k in 1:length(p)) {
    for (l in 1:length(prop_zero)) {
      for (s in 1:length(snr)){
        for (i in 1:n_iter){
          gen_data(n[j], p[k], prop_zero[l], snr[s], 
                   "~/Documents/GitHub/thesis/data/simulated/",
                   paste0("sim",'_',n[j],'_',p[k],"_",prop_zero[l],
                          "_",snr[s],"_",i))
        }
      }
    }
  }
}

