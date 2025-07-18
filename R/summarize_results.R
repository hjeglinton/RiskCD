setwd("/Users/seehanahtang/Downloads/my_research/honors_thesis")
library(tidyverse)
library(ggpubr)
library(kableExtra)


process_results <- function(df = NULL, filepath = NULL) {
  
  if (!is.null(filepath) & is.null(df)) {
    
    results <- read.csv(filepath)
    
  } else if (!is.null(df)) {
    
    results <- df
    
  }
  results <- results %>%
    separate(col = 1,
             into = c("sim", "n_data", "p_data", "prop_zero", "snr", "iter", "end"),
             sep = "_", convert = TRUE, remove = FALSE) %>%
    select(-c(sim, end, n_data, p_data))
  
  return(results)
  
}

summarize_metrics <- function(df, model) {
  
  summary <- df %>%
    filter(method == model) %>%
    group_by(n, p, prop_zero) %>%
    summarize(AUC = round(mean(auc_test), 3),
              sd = round(sd(auc_test), 3),
              brier = round(mean(brier_test), 3),
              num_nonzeros = round(mean(num_nonzeros),1),
              seconds = round(mean(seconds),2))
  
  return(summary)
  
}

percent_best <- function(results_df, methods) {
  
  results_compare <- results_df %>%
    select(n, p, prop_zero, snr, iter, method, auc_test) %>%
    filter(method %in% methods) %>%
    pivot_wider(names_from = "method", values_from = "auc_test")
  
  method_index <- seq(6, ncol(results_compare))
  
  results_compare[paste0(names(results_compare)[method_index], "_best")] <- 0
  
  for (i in 1:nrow(results_compare)) {
    
    max_auc <- max(round(results_compare[i, method_index],4), na.rm = T)
    max_method <- which(round(results_compare[i, method_index],4) == max_auc)
    
    results_compare[i, max_method + max(method_index)] = 1
    
  }
  
  return(results_compare)
  
}

####################################

# Read in results
# results <- process_results(filepath = "RiskCD/results/simulated/results_range3.csv")
results <- process_results(filepath = "RiskCD/results/simulated/results_range5.csv")
results <- results [(!(results$method=="RiskCD") & !(results$method=="RiskCD-CV")),]

results$method <- factor(results$method, levels = 
                           c("LR", "Rounded LR", "Lasso", "Rounded Lasso", "FasterRisk", "FasterRisk-CV", 
                                                    "AnnealScore", "AnnealScore-CV"))

#### Table -- no CV ####

# %Best 
percent_best_nocv_df <- percent_best(results,
                                     methods = c("FasterRisk", "AnnealScore", "Rounded LR")) %>%
  group_by(n, p, prop_zero) %>%
  summarize(AS_best = mean(AnnealScore_best)*100,
            rLR_best = mean(`Rounded LR_best`)*100,
            FR_best = mean(FasterRisk_best)*100)

percent_best_cv_df <- percent_best(results,
                                methods = c("FasterRisk-CV", "AnnealScore-CV", "Rounded Lasso")) %>%
  group_by(n, p, prop_zero) %>%
  summarize(AS_CV_best = mean(`AnnealScore-CV_best`)*100,
            rLasso_best = mean(`Rounded Lasso_best`)*100,
            FR_CV_best = mean(`FasterRisk-CV_best`)*100)

# write.csv(percent_best_cv_df, "results/percent_best_cv.csv", row.names = FALSE)
# write.csv(percent_best_nocv_df, "results/percent_best_no_cv.csv", row.names = FALSE)


#### Table -- CV ####

# Summaries -- CV
LR_summary <- summarize_metrics(results, "LR")
rLR_summary <- summarize_metrics(results, "Rounded LR")
FR_summary <- summarize_metrics(results, "FasterRisk")
AS_summary <- summarize_metrics(results, "AnnealScore")
Lasso_summary <- summarize_metrics(results, "Lasso")
rLasso_summary <- summarize_metrics(results, "Rounded Lasso")
FR_CV_summary <- summarize_metrics(results, "FasterRisk-CV")
AS_CV_summary <- summarize_metrics(results, "AnnealScore-CV")

# rLasso vs AnnealScore-CV
round((rLasso_summary$AUC - AS_CV_summary$AUC) / rLasso_summary$AUC, 3)

# rLasso vs SA-CV
round((rLasso_summary$AUC - AS_CV_summary$AUC) / rLasso_summary$AUC, 3)

# AS-CV vs FR-CV
round((AS_CV_summary$AUC - FR_CV_summary$AUC) / AS_CV_summary$AUC, 3)


#### Accuracy plot ####

for (i in 1:nrow(results)) {
  
  coef_file_name <- str_replace(results$data[i], "data", "coef")
  coef_df <- read.csv(paste0("simulated-new/", coef_file_name))[-1,]
  
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

# --- 50 % noise -------------------------------------------------------------
nonzero_50 <- results %>%                       # build the summary table
  filter(prop_zero == 0.5) %>%                 
  group_by(method) %>%                         
  summarize(mean_acc = mean(nonzero_accuracy, na.rm = TRUE),
            sd_acc   = sd(nonzero_accuracy,   na.rm = TRUE)) 

accuracy_plot1 <- ggplot(nonzero_50) +
  geom_bar(aes(x = method, y = mean_acc), stat = "identity",
           width = 0.2, fill = "white", color = "black", linewidth = 1) +
  geom_errorbar(aes(x = method,
                    ymin = mean_acc - sd_acc/2,
                    ymax = mean_acc + sd_acc/2),
                width = 0.1, linewidth = 0.75) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = "Accuracy",
       title = "Accuracy in Detecting 50% Noise Predictors") +
  theme_bw(base_size = 20)

# --- 0 % noise --------------------------------------------------------------
nonzero_0 <- results %>%
  filter(prop_zero == 0) %>% 
  group_by(method) %>% 
  summarize(mean_acc = mean(nonzero_accuracy, na.rm = TRUE),
            sd_acc   = sd(nonzero_accuracy,   na.rm = TRUE)) %>%

accuracy_plot2 <- ggplot(nonzero_0) +
  geom_bar(aes(x = method, y = mean_acc), stat = "identity",
           width = 0.2, fill = "white", color = "black", linewidth = 1) +
  geom_errorbar(aes(x = method,
                    ymin = mean_acc - sd_acc/2,
                    ymax = mean_acc + sd_acc/2),
                width = 0.1, linewidth = 0.75) +
  geom_hline(yintercept = 0) +
  labs(x = NULL, y = "Accuracy",
       title = "Accuracy in Detecting 0% Noise Predictors") +
  theme_bw(base_size = 20)

# Combine and save -----------------------------------------------------------
both_sparsity_plots <- ggarrange(accuracy_plot1, accuracy_plot2, ncol = 1)
ggsave(both_sparsity_plots,
       filename = "results/sim_nonzero_accuracy_plots.png",
       width = 15, height = 10, dpi = 300)



#### Time plots ####

time_results <- results %>%
  select(n, p, method, seconds) %>%
  filter(method %in% c("Rounded LR", "AnnealScore", "FasterRisk", "Rounded Lasso", "AnnealScore-CV", "FasterRisk-CV"))

# no CV 
time_plot_noCV_n <- time_results %>%
  filter(method %in% c("FasterRisk", "AnnealScore")) %>%
  group_by(n, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = n, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = n, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = n, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 100, alpha = 0.8) +
  labs(x = "Number of observations (n)", y = "Seconds", linetype = "", shape = "") +
  theme_bw() +
  # scale_x_continuous(breaks = c(100, 500, 1000)) +
  scale_y_continuous(limits = c(0, 120)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

time_plot_noCV_p <- time_results %>%
  filter(method %in% c("FasterRisk", "AnnealScore")) %>%
  group_by(p, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = p, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = p, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = p, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 0.8, alpha = 0.8) +
  labs(x = "Number of candidate covariates (p)", y = "", linetype = "", shape = "") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 120)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

both_plots_noCV <- ggarrange(time_plot_noCV_n, time_plot_noCV_p, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(both_plots_noCV, width = 15, height = 8, dpi = 300,
       filename = "results/sim_time_noCV.png")

# CV
time_plot_CV_p <- time_results %>%
  filter(method %in% c("FasterRisk-CV","AnnealScore-CV")) %>%
  group_by(p, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = p, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = p, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = p, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 0.8, alpha = 0.8) +
  labs(x = "Number of candidate covariates (p)", y = "", linetype = "", shape = "") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

ggsave(time_plot_CV_p, width = 8, height = 8, dpi = 300,
       filename = "results/sim_time_CV.png")

time_summary <- time_results %>%
  group_by(n, p, method) %>%
  summarize(avg_time = round(mean(seconds, na.rm = TRUE), 2), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(n, p),
    names_from = method,
    values_from = avg_time
  ) %>%
  select(
    n, p,
    `Rounded LR`, AnnealScore, FasterRisk,       # Non-regularized
    `Rounded Lasso`, `AnnealScore-CV`, `FasterRisk-CV`  # Regularized
  )


#### Summary tables ####

table_summary_no_cv <- LR_summary %>%
  rename(LR_AUC          = AUC,
         LR_Brier        = brier,
         LR_SD           = sd,
         LR_Time         = seconds) %>%
  left_join(
    rLR_summary %>%
      rename(rLR_AUC          = AUC,
             rLR_Brier        = brier,
             rLR_SD           = sd,
             rLR_Time         = seconds),
    by = c("n", "p", "prop_zero")
  ) %>%
  left_join(
    FR_summary %>%
      rename(FR_AUC           = AUC,
             FR_Brier         = brier,
             FR_SD           = sd,
             FR_Time          = seconds),
    by = c("n", "p", "prop_zero")
  ) %>%
  left_join(
    AS_summary %>%
      rename(AS_AUC           = AUC,
             AS_Brier         = brier,
             AS_SD            = sd,
             AS_Time          = seconds),
    by = c("n", "p", "prop_zero")
  ) %>%
  left_join(
    percent_best_nocv_df,
    by = c("n", "p", "prop_zero")
  )

table_summary_cv <- Lasso_summary %>%
      rename(Lasso_AUC          = AUC,
             Lasso_Brier        = brier,
             Lasso_SD           = sd,
             Lasso_Time         = seconds) %>%

  # Join Rounded Lasso summaries
  left_join(
    rLasso_summary %>%
      rename(rLasso_AUC          = AUC,
             rLasso_Brier        = brier,
             rLasso_SD           = sd,
             rLasso_Time         = seconds),
    by = c("n", "p", "prop_zero")
  ) %>%
  
  # Join FasterRisk-CV summaries
  left_join(
    FR_CV_summary %>%
      rename(FR_CV_AUC           = AUC,
             FR_CV_Brier         = brier,
             FR_CV_SD            = sd,
             FR_CV_Time          = seconds),
    by = c("n", "p", "prop_zero")
  ) %>%
  
  # Join AnnealScore-CV
  left_join(
    AS_CV_summary %>%
      rename(AS_CV_AUC           = AUC,
             AS_CV_Brier         = brier,
             AS_CV_SD            = sd,
             AS_CV_Time          = seconds),
    by = c("n", "p", "prop_zero")
  ) %>%

  left_join(
    percent_best_cv_df,
    by = c("n", "p", "prop_zero")
  )

df_cv <- table_summary_cv %>%
  select(n, p, prop_zero,
         Lasso_AUC, prop_zero,
         rLasso_AUC, rLasso_SD, rLasso_best,
         AS_CV_AUC, AS_CV_SD, AS_CV_best,
         FR_CV_AUC, FR_CV_SD, FR_CV_best)

df_no_cv <- table_summary_no_cv %>%
  select(n, p, prop_zero,
    LR_AUC,
    rLR_AUC, rLR_SD, rLR_best,
    AS_AUC, AS_SD, AS_best,
    FR_AUC, FR_SD, FR_best
  )

write.csv(df_no_cv, "results/summary_no_cv.csv", row.names = FALSE)
write.csv(df_cv, "results/summary_cv.csv", row.names = FALSE)