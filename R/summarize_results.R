setwd("/Users/hannaheglinton/Documents/GitHub/thesis/R")

results <- read.csv("../results/simulated/results_sim_SNR_0202.csv") %>%
  select(-c(n,p)) %>%
  separate(col = 1, 
           into = c("sim", "n", "p", "prop_zero", "snr", "iter", "data"),
           sep = "_", convert = TRUE) 

View(results %>%
       group_by(p, prop_zero) %>%
       summarize(mean(nonzeros)))


summary1 <- results %>%
  filter(method %in% c("LR", "Rounded LR", "FasterRisk", "NLLCD")) %>%
  group_by(n, p, prop_zero, method) %>%
  summarize(seconds_mean = mean(seconds),
            seconds_sd = sd(seconds),
            nonzeros_mean = mean(nonzeros),
            nonzeros_sd = sd(nonzeros),
            auc_test_mean = mean(auc_test),
            auc_test_sd = sd(auc_test)) %>%
  mutate(across(seconds_mean:auc_test_sd, \(x) round(x,3))) %>%
  unite("seconds", 5:6, sep = " (") %>%
  unite("nonzeros", 6:7, sep = " (") %>%
  unite("auc", 7:8, sep = " (")




integer_methods <- c("FasterRisk", "NLLCD", "Rounded Lasso", "Rounded LR", "NLLCD with CV (lambda_min)")

summary %>%
  filter(method %in% integer_methods)


perc_best <- summary %>%
  filter(method != "LR" & method != "Lasso") %>%
  arrange(desc(auc_test)) %>%
  group_by(n, p, prop_zero, snr) %>%
  slice(1) %>%
  group_by(n, p) %>%
  summarize(FasterRisk = round(100*mean(method == "FasterRisk"),0), 
            NLLCD = round(100*mean(method == "NLLCD" | 
                                   method == "NLLCD with CV (lambda_1se)" |
                                   method == "NLLCD with CV (lambda_min)"),0),
            `Rounded Lasso` = round(100*mean(method == "Rounded Lasso"),0),
            `Rounded LR` = round(100*mean(method == "Rounded LR"), 0)) %>%
  pivot_longer(3:6, names_to = "method", values_to = "perc_best")



perc_best %>%
  filter()



FRvsNLLCD <- results %>%
  filter(method == "FasterRisk" | method == "NLLCD") %>%
  arrange(desc(auc_test)) %>%
  group_by(n, p, prop_zero, snr, iter) %>%
  slice(1) %>%
  #group_by(n, p) %>%
  summarize(winner = ifelse(mean(method == "FasterRisk") == 1, "FasterRisk", "NLLCD"),
            seconds = mean(seconds))


FRvsNLLCD %>%
  group_by(winner) %>%
  summarize(avg_n = mean(as.numeric(n)), 
            avg_p = mean(as.numeric(p)), 
            avg_prop_zero = mean(as.numeric(prop_zero)), 
            avg_snr = mean(as.numeric(snr)),
            avg_seconds = mean(seconds))


View(FRvsNLLCD %>%
  group_by(n, p, prop_zero, snr) %>%
  summarize(perc_FR = mean(winner == "FasterRisk")))


summary %>%
  filter(method %in% integer_methods) %>%
  group_by(prop_zero, method) %>%
  summarize(auc_test = mean(auc_test)) %>%
ggplot() +
  geom_point(aes(x = prop_zero, y = auc_test, color = method)) 


#write.csv(summary, "../sim_data/results_summary_0930.csv")

summary %>%
ggplot() +
  geom_violin(aes(x = method, y = auc_test))




