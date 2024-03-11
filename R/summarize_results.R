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
  select(n, p, prop_zero, snr, iter, method, auc_test) %>%
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

results_compare %>%
  group_by(n, p) %>%
  summarize(nllcd_best = mean(best == 1)*100,
            rLR_best = mean(best == 2)*100,
            FR_best = mean(best == 3)*100)


prop_best_noCV <- results %>%
  filter(method %in% c("Rounded LR", "FasterRisk", "NLLCD")) %>%
  group_by(n, p, prop_zero, snr, iter) %>%
  arrange(desc(auc_test)) %>%
  slice(1) %>%
  group_by(n, p) %>%
  summarize(FR_propbest = round(100*mean(method == "FasterRisk"),0),
            NLLCD_propbest = round(100*mean(method == "NLLCD"),0),
            roundedLR_propbest = round(100*mean(method == "Rounded LR"), 0))

prop_best_CV <- results %>%
  filter(method %in% c("Rounded Lasso", "NLLCD with CV (lambda_min)")) %>%
  group_by(n, p, prop_zero, snr, iter) %>%
  arrange(desc(auc_test)) %>%
  slice(1) %>%
  group_by(n, p) %>%
  summarize(roundedLasso_propbest = round(100*mean(method == "Rounded Lasso"),0),
            NLLCD_CV_propbest = round(100*mean(method == "NLLCD with CV (lambda_min)"),0))



RoundedLR_summary <- results %>%
  filter(method == "Rounded LR") %>%
  group_by(n, p, method) %>%
  summarize(roundedLR_AUC = round(mean(auc_test), 3), roundedLR_num_nonzeros = round(mean(num_nonzeros),1)) %>%
  select(-method)

FasterRisk_summary <- results %>%
  filter(method == "FasterRisk") %>%
  group_by(n, p, method) %>%
  summarize(FR_AUC = round(mean(auc_test), 3), FR_num_nonzeros = round(mean(num_nonzeros),1)) %>%
  select(-method)

RoundedLasso_summary <- results %>%
  filter(method == "Rounded Lasso") %>%
  group_by(n, p, method) %>%
  summarize(RoundedLasso_AUC = round(mean(auc_test), 3), RoundedLasso_num_nonzeros = round(mean(num_nonzeros),1)) %>%
  select(-method)

NLLCD_summary <- results %>%
  filter(method == "NLLCD") %>%
  group_by(n, p, method) %>%
  summarize(NLLCD_AUC = round(mean(auc_test), 3), NLLCD_num_nonzeros = round(mean(num_nonzeros),1)) %>%
  select(-method)

NLLCD_CV_summary <- results %>%
  filter(method == "NLLCD with CV (lambda_min)") %>%
  group_by(n, p, method) %>%
  summarize(NLLCD_CV_AUC = round(mean(auc_test), 3), NLLCD_CV_num_nonzeros = round(mean(num_nonzeros),1)) %>%
  select(-method)

summary_noCV <- left_join(RoundedLR_summary, FasterRisk_summary) %>%
  left_join(NLLCD_summary) %>%
  left_join(prop_best_noCV) %>%
  left_join(LR_AUC) %>%
  select(n, p, LR_AUC,
         roundedLR_AUC, roundedLR_num_nonzeros, roundedLR_propbest,
         FR_AUC, FR_num_nonzeros, FR_propbest,
         NLLCD_AUC, NLLCD_num_nonzeros, NLLCD_propbest)

summary_CV <- left_join(Lasso_AUC, RoundedLasso_summary) %>%
  left_join(NLLCD_CV_summary) %>%
  left_join(prop_best_CV) %>%
  select(n, p, Lasso_AUC,
         RoundedLasso_AUC, RoundedLasso_num_nonzeros, roundedLasso_propbest,
         NLLCD_CV_AUC, NLLCD_CV_num_nonzeros, NLLCD_CV_propbest)

write.csv(summary_noCV, "../results/simulated/summary_noCV.csv")
write.csv(summary_CV, "../results/simulated/summary_CV.csv")


