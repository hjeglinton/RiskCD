setwd("/Users/hannaheglinton/Documents/GitHub/thesis/R")

results_noCV <- read.csv("../results/simulated/results_sim_0309_noCV.csv") %>%
  select(-c(n,p)) %>%
  separate(col = 1,
           into = c("sim", "n", "p", "prop_zero", "snr", "iter", "data"),
           sep = "_", convert = TRUE)

results_FR_noCV <- read.csv("../results/simulated/results_sim_0309_FR_noCV.csv") %>%
  select(-c(n,p)) %>%
  separate(col = 1,
           into = c("sim", "n", "p", "prop_zero", "snr", "iter", "data"),
           sep = "_", convert = TRUE)

results_CD_CV <- read.csv("../results/simulated/results_sim_0309_CV.csv") %>%
  select(-c(n,p)) %>%
  separate(col = 1,
           into = c("sim", "n", "p", "prop_zero", "snr", "iter", "data"),
           sep = "_", convert = TRUE)

p10 <- results_noCV %>%
  filter(p == 10) %>%
  group_by(n, snr, method) %>%
  summarize(auc = mean(auc_test)) %>%
  pivot_wider(names_from = "method", values_from = "auc")

summarize_metrics <- function(df, model) {

  summary <- df %>%
    filter(method == model) %>%
    group_by(n, p) %>%
    summarize(AUC = round(mean(auc_test), 3),
              brier = round(mean(brier_test), 3),
              num_nonzeros = round(mean(num_nonzeros),1))

  return(summary)

}

# Summaries -- no CV
LR_summary <- summarize_metrics(results_noCV, "LR")
rLR_summary <- summarize_metrics(results_noCV, "Rounded LR")
FR_summary <- summarize_metrics(results_FR_noCV, "FasterRisk")
CD_summary <- summarize_metrics(results_noCV, "NLLCD")

# FR vs RiskCD
round((FR_summary$AUC - CD_summary$AUC) / FR_summary$AUC, 3)

# %Best -- no CV
results_compare <- bind_rows(results_noCV, results_FR_noCV) %>%
  select(n, p, prop_zero, snr, iter, method, auc_test, num_nonzeros) %>%
  filter(method %in% c("FasterRisk", "NLLCD", "Rounded LR")) %>%
  pivot_wider(names_from = "method", values_from = "auc_test") %>%
  mutate(best = NA)

for (i in 1:nrow(results_compare)) {

  max_auc <- max(round(results_compare[i, 6:8],4))
  max_method <- which(round(results_compare[i, 6:8],4) == max_auc)

  if (length(max_method) > 1) {
    results_compare$best[i] = 0
  } else {
    results_compare$best[i] = max_method
  }

}

# Summaries -- CV
lasso_summary <- summarize_metrics(results_noCV, "Lasso")
rlasso_summary <- summarize_metrics(results_noCV, "Rounded Lasso")
cdcv_summary <- summarize_metrics(results_CD_CV, "NLLCD with CV (lambda_min)")

results_compare %>%
  group_by(n, p) %>%
  summarize(nllcd_best = mean(best == 1)*100,
            rLR_best = mean(best == 2)*100,
            FR_best = mean(best == 3)*100)

# rLasso vs CD-CV
round((rlasso_summary$AUC - cdcv_summary$AUC) / rlasso_summary$AUC, 3)


# %Best -- CV
results_compare_cv <- bind_rows(results_noCV, results_CD_CV) %>%
  select(n, p, prop_zero, snr, iter, method, auc_test) %>%
  filter(method %in% c("Rounded Lasso", "NLLCD with CV (lambda_min)")) %>%
  pivot_wider(names_from = "method", values_from = "auc_test") %>%
  mutate(best = NA)

for (i in 1:nrow(results_compare_cv)) {

  max_auc <- max(round(results_compare_cv[i, 6:7],4))
  max_method <- which(round(results_compare_cv[i, 6:7],4) == max_auc)

  if (length(max_method) > 1) {
    results_compare_cv$best[i] = 0
  } else {
    results_compare_cv$best[i] = max_method
  }

}

results_compare_cv %>%
  group_by(n, p) %>%
  summarize(rlasso_best = mean(best == 1)*100,
            cdcv_best = mean(best == 2)*100)







