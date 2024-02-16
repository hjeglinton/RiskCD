setwd("/Users/hannaheglinton/Documents/GitHub/thesis/R")

results <- read.csv("../results/simulated/results_sim_SNR_0202.csv") %>%
  select(-c(n,p)) %>%
  separate(col = 1, 
           into = c("sim", "n", "p", "prop_zero", "snr", "iter", "data"),
           sep = "_", convert = TRUE) 


LR_AUC <- results %>%
  filter(method == "LR") %>%
  group_by(n, p) %>%
  summarize(LR_AUC = round(mean(auc_test),3))

Lasso_AUC <- results %>%
  filter(method == "Lasso") %>%
  group_by(n, p) %>%
  summarize(Lasso_AUC = round(mean(auc_test),3))

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
  summarize(roundedLR_AUC = round(mean(auc_test), 3), roundedLR_nonzeros = round(mean(nonzeros),1)) %>%
  select(-method) 

FasterRisk_summary <- results %>%
  filter(method == "FasterRisk") %>%
  group_by(n, p, method) %>%
  summarize(FR_AUC = round(mean(auc_test), 3), FR_nonzeros = round(mean(nonzeros),1)) %>%
  select(-method) 

RoundedLasso_summary <- results %>%
  filter(method == "Rounded Lasso") %>%
  group_by(n, p, method) %>%
  summarize(RoundedLasso_AUC = round(mean(auc_test), 3), RoundedLasso_nonzeros = round(mean(nonzeros),1)) %>%
  select(-method) 

NLLCD_summary <- results %>%
  filter(method == "NLLCD") %>%
  group_by(n, p, method) %>%
  summarize(NLLCD_AUC = round(mean(auc_test), 3), NLLCD_nonzeros = round(mean(nonzeros),1)) %>%
  select(-method) 

NLLCD_CV_summary <- results %>%
  filter(method == "NLLCD with CV (lambda_min)") %>%
  group_by(n, p, method) %>%
  summarize(NLLCD_CV_AUC = round(mean(auc_test), 3), NLLCD_CV_nonzeros = round(mean(nonzeros),1)) %>%
  select(-method) 

summary_noCV <- left_join(RoundedLR_summary, FasterRisk_summary) %>%
  left_join(NLLCD_summary) %>%
  left_join(prop_best_noCV) %>%
  left_join(LR_AUC) %>%
  select(n, p, LR_AUC, 
         roundedLR_AUC, roundedLR_nonzeros, roundedLR_propbest,
         FR_AUC, FR_nonzeros, FR_propbest,
         NLLCD_AUC, NLLCD_nonzeros, NLLCD_propbest)

summary_CV <- left_join(Lasso_AUC, RoundedLasso_summary) %>%
  left_join(NLLCD_CV_summary) %>%
  left_join(prop_best_CV) %>%
  select(n, p, Lasso_AUC, 
         RoundedLasso_AUC, RoundedLasso_nonzeros, roundedLasso_propbest,
         NLLCD_CV_AUC, NLLCD_CV_nonzeros, NLLCD_CV_propbest)

write.csv(summary_noCV, "../results/simulated/summary_noCV.csv")
write.csv(summary_CV, "../results/simulated/summary_CV.csv")


