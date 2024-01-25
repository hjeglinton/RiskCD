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
sim_data <- function(n, p, prop_zero, eps = 0){
  
  # Simulate covariates
  x <- matrix(0, nrow = n, ncol = p)
  for (i in 1:p){
    x[,i] <- rbinom(n, 1, runif(1, 0.1, 0.9))
  }
  
  # Simulate coefficients
  p_zero <- p*floor(prop_zero)
  p_nonzero <- p - p_zero
  
  coef_small <- runif(floor(p_nonzero*(3/4)), 0.5, 1.5)
  coef_large <- runif(ceiling(p_nonzero*(1/4)), 3, 5)
  coef <- c(coef_small, coef_large, rep(0, p_zero))
  
  # Get intercept
  vals <- x %*% coef
  intercept <- -mean(vals)
  
  # Simulate outcome
  vals <- vals + intercept + rnorm(n, 0, eps * sd(vals))
  probs <- exp(vals)/(1 + exp(vals))
  y <- rbinom(n, 1, prob = probs)
  
  return(list(x=x,y=y,coef=coef, intercept =intercept))
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
gen_data <- function(n, p, prop_zero, eps, folder, filename){

  # Simulate data
  data <- sim_data(n, p, prop_zero, eps)
  
  # data file
  df <- as.data.frame(data$x)
  df$y <- data$y
  df <- df %>%
    select(y, everything())
  write.csv(df, paste0(folder, filename,"_data.csv"), row.names=FALSE)
  
  # coefficients file
  coef_df <- data.frame(names = c("Intercept", names(df)[2:length(df)]),
                        vals = c(data$intercept, data$coef))
  write.csv(coef_df, paste0(folder,filename,"_coef.csv"), row.names=FALSE)
  
}

# Example of generating data and saving
set.seed(5)
n = c(100, 500, 1000, 5000)
p = c(10, 25, 50)
prop_zero = c(0, 0.25, 0.50, 0.75)
eps = c(0, 0.2, 0.4, 0.6, 0.8)

for (j in 1:length(n)){
  for (k in 1:length(p)) {
    for (l in 1:length(prop_zero)) {
      for (e in 1:length(eps)){
        for (i in 1:5){
          gen_data(n[j], p[k], prop_zero[l], eps[e], 
                   "~/Documents/GitHub/thesis/data/simulated/",
                   paste0("sim",'_',n[j],'_',p[k],"_",prop_zero[l],
                          "_",as.integer(10*eps[e]),"_",i))
        }
      }
    }
  }
}

