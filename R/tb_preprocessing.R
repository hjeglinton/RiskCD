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

mode_imputation <- function(x) {
  ux <- unique(x)
  mode <- ux[which.max(tabulate(match(x, ux)))]
  
  new_x <- replace_na(x, mode)
  return(new_x)
}

kharitode_preprocessing <- function(tb_df) {
  # Select variables and rename
  tb_df <- tb_df %>% select(c(xpert_status_fac, age_group, sex, hiv_status_fac,
                              other_conditions_fac___3, symp_fac___1, symp_fac___2,
                              symp_fac___3, symp_fac___4, length_symp_wk_fac, 
                              length_symp_days_fac, length_symp_unit_fac, smk_fac, 
                              dx_tb_past_fac, educ_fac)) %>%
    rename(tb = xpert_status_fac, hiv_pos = hiv_status_fac,
           diabetes = other_conditions_fac___3, cough = symp_fac___1,
           fever = symp_fac___2, weight_loss = symp_fac___3, 
           night_sweats = symp_fac___4, ever_smoke = smk_fac, past_tb = dx_tb_past_fac,
           education = educ_fac)
  
  tb_df$tb <- case_when(tb_df$tb == 1 ~ 1, tb_df$tb == 2 ~ 0)
  tb_df$male <- case_when(tb_df$sex == 1 ~ 1, tb_df$sex == 2 ~ 0)
  tb_df$hiv_pos <- case_when(tb_df$hiv_pos == 1 ~ 1, tb_df$hiv_pos %in% c(2,77) ~ 0)
  tb_df$ever_smoke <- case_when(tb_df$ever_smoke %in% c(1,2) ~ 1, 
                                tb_df$ever_smoke == 3 ~ 0)
  tb_df$past_tb <- case_when(tb_df$past_tb == 1 ~ 1, tb_df$past_tb == 2 ~ 0)
  
  # Categorize education to HS and above
  tb_df$hs_less <- case_when(tb_df$education <= 12 ~ 1,
                             tb_df$education %in% c(13, 14) ~ 0,
                             TRUE ~ NA)
  tb_df <- tb_df %>% select(-c(education))
  
  # Categorize number of weeks
  tb_df$two_weeks_symp <- case_when(tb_df$length_symp_days_fac <= 14 | 
                                      tb_df$length_symp_wk_fac < 2 ~ 0, 
                                    tb_df$length_symp_unit_fac == 77 |
                                      tb_df$length_symp_wk_fac == 999 ~ NA,
                                    TRUE ~ 1)
  tb_df <- tb_df %>% select(-c(length_symp_wk_fac, length_symp_days_fac,
                               length_symp_unit_fac, sex))
  
  # Total number of symptoms
  tb_df$num_symptoms <- tb_df$fever + tb_df$weight_loss +
    tb_df$cough+tb_df$night_sweats
  tb_df <- tb_df %>% select(-c(night_sweats, weight_loss, cough, fever))
  
  # Exclude participants with no TB symptoms
  tb_df <- tb_df %>%
    filter(num_symptoms != 0)
  
  # All factors
  tb_df[] <- lapply(tb_df, function(x){return(as.factor(x))})
  tb_df$age_group <- relevel(tb_df$age_group, ref="[55,99)")

  # Single imputation (mode)
  tb_df <- tb_df %>%
    mutate_all(mode_imputation)
  
  return(tb_df)
  
}

stomp_preprocessing <- function(df) {
  
  # Select variables and rename
  df <- df %>%
    select(why_a_case_xpert, com_why_a_case_xpert,
           why_a_ctrl_xpert, com_why_a_ctrl_xpert,
           age_group, sex_female, hivcat1, medical_conditions___3,
           cough = symptoms_past_week___1, fever = symptoms_past_week___3,
           night_sweats = symptoms_past_week___5, weight_loss, 
           smoked_tobacco_in_past, treated_past, highest_grade_education,
           cough_weeks, weightloss_weeks, nightsweats_weeks, feverchills_weeks) %>%
    
    rename(hiv_pos = hivcat1, diabetes = medical_conditions___3, 
           ever_smoke = smoked_tobacco_in_past, past_tb = treated_past) %>%
    
    mutate(tb = case_when(why_a_case_xpert == 1 | com_why_a_case_xpert == 1 ~ 1,
                          why_a_ctrl_xpert == 2 | com_why_a_ctrl_xpert == 2 ~ 0),
           male = case_when(sex_female == 1 ~ 0,
                            sex_female == 0 ~ 1),
           hs_less = case_when(highest_grade_education < 5 ~ 1,
                               highest_grade_education >= 5 ~ 0),
           num_symptoms = cough + fever + night_sweats + weight_loss,
           two_weeks_symp = case_when(cough_weeks >= 2 | weightloss_weeks >= 2 |
                             nightsweats_weeks >= 2 | feverchills_weeks >= 2 ~ 1,
                             TRUE ~ 0)) %>%
    
    select(tb, age_group, hiv_pos, diabetes, ever_smoke, past_tb, male, hs_less,
           two_weeks_symp, num_symptoms)
  
  
  # All factors
  df[] <- lapply(df, function(x){return(as.factor(x))})
  df$age_group <- relevel(df$age_group, ref="[55,99)")
  
  # Single imputation (mode)
  df <- df %>%
    mutate_all(mode_imputation)
  
  return(df)
  
  
}


