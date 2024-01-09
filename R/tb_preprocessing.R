suppressMessages(library(tidyverse))

tb_preprocessing <- function(tb_df) {
  
  # we choose to look at whether 100% of doses are taken on time 
  # this is a fairly balanced outcome
  tb_df$adherence_outcome <- (tb_df$PCTadherence_sensi < 100)
  
  # drop X, patient id, and other variables we don't want to include in our model
  tb_df <- tb_df %>% select(-c(X, # index 
                               prop_misseddoses_int, # not in dictionary
                               PCTadherence, 
                               PCTadherence_sensi, 
                               PTID2, # study ID
                               hunger_freq, # "not including in the model"
                               health_ctr, # 71 health centers 
                               post_tb_pt_work)) #  not in dictionary 
  
  # family support - simplify by taking median response to get almost integer value
  fam_vars <- c("fam_affection", "fam_parent_getalong", "fam_others_getalong", 
                "fam_confide", "fam_parent_emosupport", "fam_others_emosupport", 
                "fam_support_others", "fam_satisfied")
  tb_df$family_median <- apply(tb_df[, fam_vars], 1, median, na.rm=TRUE)
  tb_df <- tb_df %>% select(-all_of(c(fam_vars, "fam_support")))
  
  # evaluation of health services
  health_serv_vars <- c("aten_wait", "aten_respect", "aten_explain", "aten_space",
                        "aten_concern","aten_satis_hours")
  tb_df$health_svc_median <- apply(tb_df[, health_serv_vars], 
                                   1, median, na.rm=TRUE)
  tb_df <- tb_df %>% select(-all_of(c(health_serv_vars, "healthsvc_satis")))
  
  # motivation
  motivation_vars <- c("motiv_no_transmit", "motiv_fam_worry", "motiv_study_work",
                       "motiv_activities")
  tb_df$motiv_median <- apply(tb_df[, motivation_vars], 1, median, na.rm=TRUE)
  tb_df <- tb_df %>% select(-all_of(c(motivation_vars,  "motiv_summary")))
  
  # tb disinformation
  knowledge_vars <- c("conoc_cure", "conoc_missed_doses", "conoc_default")   
  tb_df$knowledge_median <- apply(tb_df[, knowledge_vars], 1, median, na.rm=TRUE)
  tb_df <- tb_df %>% select(-all_of(c(knowledge_vars,  "tb_knowledge")))
  
  # continuous variables
  cont_var <- c("self_eff", "tx_mos", "audit_tot", "stig_tot", 
                "phq9_tot", "age_BLchart", "ace_score")
  
  # categorical but treat as continuous with small range and int values
  cont_cat_vars <- c("tobacco_freq", "covid_es", "pills", "adr_freq", 
                     "fam_accompany_dot", "fam_dislikefriends", 
                     "autonomy_obedient", "stig_health_ctr", "family_median",
                     "health_svc_median", "motiv_median", "knowledge_median")
  
  # pills and adr_freq need to multiplied by 4 - does not match data dictionary
  tb_df$pills <- 4*tb_df$pills
  tb_df$adr_freq <- 4*tb_df$adr_freq
  
  # change to categorical
  tb_df$current_sx_none <- as.factor(tb_df$current_sx_none) 
  tb_df$tto_anterior_tb <- as.factor(tb_df$tto_anterior_tb) 
  
  # make continuous var updates
  tb_df <- tb_df %>%
    mutate(age_cat = case_when(age_BLchart < 16 ~ "< 16", 
                               age_BLchart < 18 ~ "16-17", TRUE ~ "18+"),
           audit_cat = case_when(audit_tot==0 ~ "0", TRUE ~ ">0"),
           ace_cat = case_when(ace_score == 0 ~ "0",
                               ace_score == 1 ~ "1",
                               ace_score > 1 ~ "> 1"),
           tx_mos_cat = case_when(tx_mos <= 6 ~ "<= 6 mos", 
                                  TRUE ~ "> 6 mos")) %>%
    select(-c(self_eff, tx_mos, audit_tot, stig_tot, phq9_tot, age_BLchart,
              ace_score)) %>%
    na.omit()
  
  # some categorical variable have levels with few observations
  # drop ram and regular_drug since only 5 observations in class 1
  tb_df <- tb_df %>% select(-c(ram, regular_drug))
  
  # education levels and monitor1 could be combined but we drop since not 
  # included in data documentation and unclear best way to relevel
  tb_df <- tb_df %>% select(-c(edu_level_mom, edu_level_dad, monitor1))
  
  return(tb_df)
}

tb_as_matrix <- function(tb_df) {
  
  # model matrix and save
  tb_matrix <- model.matrix(adherence_outcome ~ ., data=tb_df)
  tb_matrix <- cbind(tb_matrix, adherence_outcome = tb_df$adherence_outcome)
  
  return(tb_matrix)
  
}

