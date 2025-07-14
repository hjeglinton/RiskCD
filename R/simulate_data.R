library(tidyverse)


#' Simulate data and coefficients
#'
#' @param n Number of observations for train set.
#' @param p Total number of variables.
#' @param p_zero Proportion of variables with a coefficient of zero.
#' @param eps Value between 0 and 1 representing amount of noise.
#' @return List with simulated covariates (x), outcome (y), coefficients (coef),
#' and intercept.
sim_data <- function(n, p, prop_zero, snr){

  # Simulate coefficients
  p_zero <- floor(p*prop_zero)
  p_nonzero <- p - p_zero
  coef <- c(runif(p_nonzero, 0, 10), rep(0, p_zero))
  
  # Total n
  n_df <- 2*n

  # Simulate covariates
  x <- matrix(0, nrow = n_df, ncol = p)
  for (i in 1:p){
    x[,i] <- rbinom(n_df, 1, runif(1, 0.1, 0.9))
  }

  # Get intercept
  vals <- x %*% coef
  intercept <- -mean(vals)

  # Get epsilon values
  sd_epsilon <- sqrt(var(vals, na.rm = TRUE) / snr)
  eps <- matrix(rnorm(n_df, mean = 0, sd = sd_epsilon), ncol = 1)

  # Simulate outcome
  vals <- vals + intercept + eps
  probs <- exp(vals)/(1 + exp(vals))
  y <- rbinom(n_df, 1, prob = probs)

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
  df$train <- c(rep(1, n), rep(0, n))
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
n = c(200)
p_prop = c(0.05, 0.1)
prop_zero = c(0.0)
snr = c(3)
n_iter = 40

for (j in 1:length(n)){
  for (k in 1:length(p_prop)) {
    for (l in 1:length(prop_zero)) {
      for (s in 1:length(snr)){
        for (i in 1:n_iter){
          p = p_prop[k]*n[j]
          gen_data(n[j], p, prop_zero[l], snr[s],
                   "/Users/alice/Dropbox/RiskCD/data/simulated-test/",
                   paste0("sim",'_',n[j],'_',p,"_",prop_zero[l],
                          "_",snr[s],"_",i))
        }
      }
    }
  }
}

