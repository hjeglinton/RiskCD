setwd("/Users/hannaheglinton/Documents/GitHub/thesis/R")
library(tidyverse)
library(ggpubr)


process_results <- function(df = NULL, filepath = NULL) {

  if (!is.null(filepath) & is.null(df)) {

    results <- read.csv(filepath)

  } else if (!is.null(df)) {

    results <- df

  }
  results <- results %>%
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
results_noCV <- process_results(filepath = "../results/simulated/results_sim_0309_noCV.csv") %>%
  filter(method %in% c("NLLCD", "Lasso", "LR"))
results_CD_CV <- process_results(filepath = "../results/simulated/results_sim_0317_CD_CV.csv")
results_FR_noCV <- process_results(filepath = "../results/simulated/results_sim_0318_FR_noCV.csv")
results_FR_CV_p10 <- process_results(filepath = "../results/simulated/results_sim_0314_FR_CV_p10.csv")
results_FR_CV_p25 <- process_results(filepath = "../results/simulated/results_sim_0314_FR_CV_p25.csv")
results_FR_CV_p50_n100 <- process_results(filepath = "../results/simulated/results_sim_0319_FR_CV_p50.csv")
results_rounded <- process_results(filepath = "../results/simulated/results_sim_0405_rounded.csv")

#### Table -- no CV ####

# Summaries -- no CV
LR_summary <- summarize_metrics(results_noCV, "LR")
rLR_summary <- summarize_metrics(results_rounded, "Rounded LR")
FR_summary <- summarize_metrics(results_FR_noCV, "FasterRisk")
CD_summary <- summarize_metrics(results_noCV, "NLLCD")

# rLR vs RiskCD
round((rLR_summary$AUC - CD_summary$AUC) / rLR_summary$AUC, 3)


# %Best -- no CV
percent_best(bind_rows(results_noCV, results_FR_noCV, results_rounded),
             methods = c("FasterRisk", "NLLCD", "Rounded LR")) %>%
  group_by(n, p) %>%
  summarize(nllcd_best = mean(NLLCD_best)*100,
            rLR_best = mean(`Rounded LR_best`)*100,
            FR_best = mean(FasterRisk_best)*100)

#### Table -- CV ####

# Summaries -- CV
lasso_summary <- summarize_metrics(results_noCV, "Lasso")
rlasso_summary <- summarize_metrics(results_rounded, "Rounded Lasso")
cdcv_summary <- summarize_metrics(results_CD_CV, "NLLCD with CV (lambda_min)")
frcv_summary <- summarize_metrics(bind_rows(results_FR_CV_p10, results_FR_CV_p25, results_FR_CV_p50_n100),
                                  "FasterRisk-CV")



# rLasso vs CD-CV
round((rlasso_summary$AUC - cdcv_summary$AUC) / rlasso_summary$AUC, 3)


# %Best -- CV
percent_best(filter(bind_rows(results_rounded, results_CD_CV, results_FR_CV_p10, results_FR_CV_p25, results_FR_CV_p50_n100)),
             methods = c("Rounded Lasso", "NLLCD with CV (lambda_min)",
                         "FasterRisk-CV")) %>%
  group_by(n, p) %>%
  summarize(rLasso_best = mean(`Rounded Lasso_best`)*100,
            riskcd_best = mean(`NLLCD with CV (lambda_min)_best`)*100,
            FR_best = mean(`FasterRisk-CV_best`)*100)


#### Accuracy plot ####

results <- bind_rows(results_CD_CV, results_FR_CV_p10,
                              results_FR_CV_p25, results_FR_noCV,
                              results_noCV, results_rounded) %>%
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
  geom_bar(aes(x = method, y = mean_acc), stat = "identity",
           width = 0.2, fill = "white", color = "black", linewidth = 1) +
  geom_errorbar(aes(x = method, ymin = mean_acc - sd_acc/2, ymax = mean_acc + sd_acc/2),
                width = 0.1, linewidth = 0.75) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Accuracy", title = "Accuracy in Detecting 50% Noise Predictors") +
  theme_bw() +
  theme(text = element_text(size = 20))

nonzero_0 <- results %>%
  filter(prop_zero == 0) %>%
  group_by(method) %>%
  summarize(mean_acc = mean(nonzero_accuracy), sd_acc = sd(nonzero_accuracy))

accuracy_plot2 <-  ggplot(nonzero_0) +
  geom_bar(aes(x = method, y = mean_acc), stat = "identity",
           width = 0.2, fill = "white", color = "black", linewidth = 1) +
  geom_errorbar(aes(x = method, ymin = mean_acc - sd_acc/2, ymax = mean_acc + sd_acc/2),
                width = 0.1, linewidth = 0.75) +
  geom_hline(yintercept = 0) +
  labs(x = "", y = "Accuracy ", title = "Accuracy in Detecting 0% Noise Predictors") +
  theme_bw() +
  theme(text = element_text(size = 20))

both_sparsity_plots <- ggarrange(accuracy_plot1, accuracy_plot2, ncol = 1)
ggsave(both_sparsity_plots, width = 15, height = 10, dpi = 300,
       filename = "../results/simulated/images/sim_nonzero_accuracy_plots.png")

#### Time plots ####


time_results <- results %>%
  select(p, n, method, seconds) %>%
  filter(method %in% c("FasterRisk", "FasterRisk-CV", "RiskCD", "RiskCD-CV"))

# no CV (includes p = 50)
time_plot_noCV_n <- time_results %>%
  filter(method %in% c("FasterRisk", "RiskCD")) %>%
  group_by(n, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = n, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = n, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = n, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 100, alpha = 0.8) +
  labs(x = "Number of observations (n)", y = "Seconds", linetype = "", shape = "") +
  theme_bw() +
  scale_x_continuous(breaks = c(100, 500, 1000, 5000)) +
  scale_y_continuous(breaks = seq(0, 400, 100), limits = c(0, 450)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

time_plot_noCV_p <- time_results %>%
  filter(method %in% c("FasterRisk", "RiskCD")) %>%
  group_by(p, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = p, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = p, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = p, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 0.8, alpha = 0.8) +
  labs(x = "Number of candidate covariates (p)", y = "", linetype = "", shape = "") +
  theme_bw() +
  scale_x_continuous(breaks = c(10, 25, 50)) +
  scale_y_continuous(breaks = seq(0, 400, 100), limits = c(0, 450)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

both_plots_noCV <- ggarrange(time_plot_noCV_n, time_plot_noCV_p, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(both_plots_noCV, width = 15, height = 8, dpi = 300,
       filename = "../results/simulated/images/sim_time_noCV.png")

# CV (excludes p = 50)
time_plot_CV_n <- time_results %>%
  filter(method %in% c("FasterRisk-CV", "RiskCD-CV", p != 50)) %>%
  group_by(n, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = n, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = n, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = n, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 100, alpha = 0.8) +
  labs(x = "Number of observations (n)", y = "Seconds", linetype = "", shape = "") +
  theme_bw() +
  scale_x_continuous(breaks = c(100, 500, 1000, 5000)) +
  scale_y_continuous(breaks = seq(0, 100, 50), limits = c(0, 100)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

time_plot_noCV_p <- time_results %>%
  filter(method %in% c("FasterRisk-CV", "RiskCD-CV")) %>%
  group_by(p, method) %>%
  summarize(avg_time = mean(seconds), sd = sd(seconds)) %>%
  ggplot() +
  geom_point(aes(x = p, y = avg_time, shape = method), size = 3) +
  geom_line(aes(x = p, y = avg_time, linetype = method)) +
  geom_errorbar(aes(x = p, ymin = avg_time - sd/2, ymax = avg_time + sd/2, linetype = method), width = 0.8, alpha = 0.8) +
  labs(x = "Number of candidate covariates (p)", y = "", linetype = "", shape = "") +
  theme_bw() +
  scale_x_continuous(breaks = c(10, 25, 50)) +
  scale_y_continuous(breaks = seq(0, 100, 50), limits = c(0, 100)) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "bottom", text = element_text(size = 20))

both_plots_noCV <- ggarrange(time_plot_noCV_n, time_plot_noCV_p, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave(both_plots_noCV, width = 15, height = 8, dpi = 300,
       filename = "../results/simulated/images/sim_time_noCV.png")
