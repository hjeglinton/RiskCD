setwd("/Users/hannaheglinton/Documents/GitHub/thesis/R")
library(tidyverse)




process_results <- function(file_path) {

  results <- read.csv(file_path) %>%
    select(-c(n,p)) %>%
    separate(col = 1,
             into = c("sim", "n", "p", "prop_zero", "snr", "iter", "end"),
             sep = "_", convert = TRUE, remove = FALSE) %>%
    select(-c(sim, end))

  return(results)

}

summarize_metrics <- function(df, model) {

  summary <- df %>%
    filter(method == model) %>%
    group_by(n, p) %>%
    summarize(AUC = round(mean(auc_test), 3),
              brier = round(mean(brier_test), 3),
              num_nonzeros = round(mean(num_nonzeros),1))

  return(summary)

}

percent_best <- function(results_df, methods) {

  results_compare <- results_df %>%
    select(n, p, prop_zero, snr, iter, method, auc_test) %>%
    filter(method %in% methods) %>%
    pivot_wider(names_from = "method", values_from = "auc_test") %>%
    mutate(best = NA)

  method_index <- seq(6, ncol(results_compare) - 1)


  for (i in 1:nrow(results_compare)) {

    max_auc <- max(round(results_compare[i, method_index],4))
    max_method <- which(round(results_compare[i, method_index],4) == max_auc)

    if (length(max_method) > 1) {
      results_compare$best[i] = 0
    } else {
      results_compare$best[i] = max_method
    }

  }

  return(results_compare)

}

# Read in results
results_noCV <- process_results("../results/simulated/results_sim_0309_noCV.csv")
results_CD_CV <- process_results("../results/simulated/results_sim_0317_CD_CV.csv")
results_FR_noCV <- process_results("../results/simulated/results_sim_0314_FR_noCV.csv")
results_FR_CV_p10 <- process_results("../results/simulated/results_sim_0314_FR_CV_p10.csv")
results_FR_CV_p25 <- process_results("../results/simulated/results_sim_0314_FR_CV_p25.csv")

#### Table -- no CV ####

# Summaries -- no CV
LR_summary <- summarize_metrics(results_noCV, "LR")
rLR_summary <- summarize_metrics(results_noCV, "Rounded LR")
FR_summary <- summarize_metrics(results_FR_noCV, "FasterRisk")
CD_summary <- summarize_metrics(results_noCV, "NLLCD")

# FR vs RiskCD
round((FR_summary$AUC - CD_summary$AUC) / FR_summary$AUC, 3)


# %Best -- no CV
percent_best(bind_rows(results_noCV, results_FR_noCV),
             methods = c("FasterRisk", "NLLCD", "Rounded LR")) %>%
  group_by(n, p) %>%
  summarize(nllcd_best = mean(best == 1)*100,
            rLR_best = mean(best == 2)*100,
            FR_best = mean(best == 3)*100)

#### Table -- CV ####

# Summaries -- CV
lasso_summary <- summarize_metrics(results_noCV, "Lasso")
rlasso_summary <- summarize_metrics(results_noCV, "Rounded Lasso")
cdcv_summary <- summarize_metrics(results_CD_CV, "NLLCD with CV (lambda_min)")
frcv_summary <- summarize_metrics(bind_rows(results_FR_CV_p10, results_FR_CV_p25), 
                                  "FasterRisk-CV")



# rLasso vs CD-CV
round((rlasso_summary$AUC - cdcv_summary$AUC) / rlasso_summary$AUC, 3)


# %Best -- CV
percent_best(filter(bind_rows(results_noCV, results_CD_CV, results_FR_CV_p25),
                    p == 25),
             methods = c("Rounded Lasso", "NLLCD with CV (lambda_min)",
                         "FasterRisk-CV")) %>%
  group_by(n, p) %>%
  summarize(rlasso_best = mean(best == 1)*100,
            cdcv_best = mean(best == 2)*100,
            frcv_best = mean(best == 3)*100)


#### Accuracy plot ####

results <- bind_rows(results_CD_CV, results_FR_CV_p10,
                              results_FR_CV_p25, results_FR_noCV,
                              results_noCV) %>%
  mutate(nonzero_accuracy = NA)

results$method[results$method == "NLLCD with CV (lambda_min)"] <- "RiskCD-CV"
results$method[results$method == "NLLCD"] <- "RiskCD"
results$method <- factor(results$method, levels = c("LR", "Rounded LR", "Lasso", "Rounded Lasso", "FasterRisk", "FasterRisk-CV", "RiskCD", "RiskCD-CV"))

for (i in 1:nrow(results)) {
  
  coef_file_name <- str_replace(results$data[i], "data", "coef")
  coef_df <- read.csv(paste0("../data/simulated/", coef_file_name))[-1,]
  
  nonzero_true <- case_when(coef_df$vals == 0 ~ 0,
                            TRUE ~ 1)
  
  observed <- str_split(results$which_nonzeros[i], "_")
  
  nonzero_obs <- seq(1, length(nonzero_true))
  nonzero_obs <- case_when(nonzero_obs %in% observed[[1]] ~ 1,
                           TRUE ~ 0)
  
  true_pos <- sum(nonzero_obs == 1 & nonzero_true == 1)
  true_neg <- sum(nonzero_obs == 0 & nonzero_true == 0)
  
  results$nonzero_accuracy[i] <- (true_pos + true_neg) / length(nonzero_true)
  
  
}

nonzero_50 <- results %>%
  filter(prop_zero == 0.5) %>%
  group_by(method) %>%
  summarize(mean_acc = mean(nonzero_accuracy, na.rm = TRUE), 
            sd_acc = sd(nonzero_accuracy, na.rm = TRUE))


accuracy_plot1 <-  ggplot(nonzero_50) + 
  geom_bar(aes(x = method, y = mean_acc), stat = "identity", width = 0.2) + 
  geom_errorbar(aes(x = method, ymin = mean_acc - sd_acc/2, ymax = mean_acc + sd_acc/2), width = 0.1) + 
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Accuracy", title = "Accuracy in Detecting 50% Noise Predictors") + 
  theme_bw() + 
  theme(text = element_text(size = 20))

nonzero_0 <- results %>%
  filter(prop_zero == 0) %>%
  group_by(method) %>%
  summarize(mean_acc = mean(nonzero_accuracy), sd_acc = sd(nonzero_accuracy))

accuracy_plot2 <-  ggplot(nonzero_0) + 
  geom_bar(aes(x = method, y = mean_acc), stat = "identity", width = 0.2) + 
  geom_errorbar(aes(x = method, ymin = mean_acc - sd_acc/2, ymax = mean_acc + sd_acc/2), width = 0.1) + 
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Accuracy ", title = "Accuracy in Detecting 0% Noise Predictors") + 
  theme_bw() + 
  theme(text = element_text(size = 20))

both_sparsity_plots <- ggarrange(accuracy_plot1, accuracy_plot2, ncol = 1)
ggsave(both_sparsity_plots, width = 15, height = 10, dpi = 300, 
       filename = "../results/simulated/sim_nonzero_accuracy_plots.png")

#### Time plot ####

time_results <- results %>%
  select(p, n, method, seconds) %>%
  filter(method %in% c("FasterRisk", "FasterRisk-CV", "RiskCD", "RiskCD-CV")) 


time_plot_n <- time_results %>%
  group_by(n, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() + 
  geom_point(aes(x = n, y = avg_time, color = method)) + 
  geom_line(aes(x = n, y = avg_time, color = method)) + 
  geom_errorbar(aes(x = n, ymin = avg_time - sd/2, ymax = avg_time + sd/2, color = method), width = 100, alpha = 0.8) + 
  labs(x = "Number of observations (n)", y = "Seconds", color = "") + 
  theme_bw() + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(100, 500, 1000, 5000)) + 
  scale_y_continuous(breaks = seq(0, 500, 100), limits = c(0, 500)) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

time_plot_n

time_plot_p <- time_results %>%
  group_by(p, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() + 
  geom_point(aes(x = p, y = avg_time, color = method)) + 
  geom_line(aes(x = p, y = avg_time, color = method)) + 
  geom_errorbar(aes(x = p, ymin = avg_time - sd/2, ymax = avg_time + sd/2, color = method), width = 0.8, alpha = 0.8) + 
  labs(x = "Number of candidate covariates (p)", y = "", color = "") + 
  theme_bw() + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = c(10, 25, 50)) + 
  scale_y_continuous(breaks = seq(0, 500, 100), limits = c(0, 500)) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

time_plot_p
