################################################################################

# Precarious employment and health - Understanding Society
# 3-02	Complete case outcome analysis (by age) – propensity score matched data 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) complete case descriptive and outcome analysis by age for PS matched data
# (b) complete case descriptive and outcome analysis by age for IPTW data


################################################################################

## remove any existing objects from global environment
rm(list=ls()) 

## turn off scientific notation
options(scipen = 999)

################################################################################
#####                            install packages                          #####
################################################################################

### use tableone v0.12.0 and survey v4.0 so that svy table works
#remotes::install_version("tableone", version = "0.12.0")
#remotes::install_version("survey", version = "4.0")

library(tidyverse) # all kinds of stuff 
library(fastDummies) # for creating dummy variables
library(Matching) # PS matching
library(MatchIt) # IPTW
library(tableone) # for creating table one
library(survey) # for PS weighting
library(reshape2) # Reorganizing data
library(broom) # for tidying regression outputs into df format
library(broom.mixed) # for tidying MLM regression outputs into df format
library(glmmTMB) # for multi-level modelling (faster than lme4)
#library(doParallel)
#library(foreach)
#cluster = makeCluster(detectCores() - 1) # convention to leave 1 core for OS
#registerDoParallel(cluster)

################################################################################
#####                         load and prepare data                        #####
################################################################################

####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")

## numeric non-normal vars
nonnorm_vec <- (c("age_dv_t0", "sf12mcs_dv_t0", "sf12pcs_dv_t0"))

####load eligible cases --------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")


### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
pair_cc_analytic$srh_bin_t0 <- factor(pair_cc_analytic$srh_bin_t0,
                                      levels = c("excellent/very good", 
                                                 "good/fair/poor"))
pair_cc_analytic$srh_bin_t1 <- factor(pair_cc_analytic$srh_bin_t1,
                                      levels = c("excellent/very good", 
                                                 "good/fair/poor"))

pair_cc_analytic$ghq_case4_t0 <- factor(pair_cc_analytic$ghq_case4_t0,
                                        levels = c("0-3", "4 or more"))
pair_cc_analytic$ghq_case4_t1 <- factor(pair_cc_analytic$ghq_case4_t1,
                                        levels = c("0-3", "4 or more"))



### convert SF-12 outcomes to numeric to allow svyglm to work
pair_cc_analytic$sf12pcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t1)
pair_cc_analytic$sf12mcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t1)



#### keep only variables required for propensity score (and other analysis) ----
pair_cc_analytic <- pair_cc_analytic %>% 
  dplyr::select(pidp, all_of(c(cov_vector, cov_vector2, outcome_vector))) #%>% 
#  dplyr::select(-c(psu, strata, wt_name, wt_value))

#### convert pidp to factor 
pair_cc_analytic$pidp <- factor(pair_cc_analytic$pidp)

#### convert scale vars to numeric ---------------------------------------------
pair_cc_analytic <- pair_cc_analytic %>% 
  #  mutate(fimnnet_dv_t0 = as.numeric(fimnnet_dv_t0)) %>%
  ## these to character first so score converted rather than factor level
  mutate(sf12mcs_dv_t0 = as.character(sf12mcs_dv_t0),
         sf12pcs_dv_t0 = as.character(sf12pcs_dv_t0),
         sf12mcs_dv_t1 = as.character(sf12mcs_dv_t1),
         sf12pcs_dv_t1 = as.character(sf12pcs_dv_t1)) %>% 
  mutate(sf12mcs_dv_t0 = as.numeric(sf12mcs_dv_t0),
         sf12pcs_dv_t0 = as.numeric(sf12pcs_dv_t0),
         sf12mcs_dv_t1 = as.numeric(sf12mcs_dv_t1),
         sf12pcs_dv_t1 = as.numeric(sf12pcs_dv_t1))

## change any character vars to factor
pair_cc_analytic[sapply(pair_cc_analytic, is.character)] <- lapply(pair_cc_analytic[sapply(pair_cc_analytic, is.character)], as.factor)

#### revel exposure vars so that unexposed is lower ref category ---------------
pair_cc_analytic$exposure1 <- factor(pair_cc_analytic$exposure1,
                                     levels = rev(levels(pair_cc_analytic$exposure1)))

pair_cc_analytic$exposure2 <- factor(pair_cc_analytic$exposure2,
                                     levels = rev(levels(pair_cc_analytic$exposure2)))

#### create dummy variables for categorical vars with>2 cats -------------------

pair_cc_analytic <- pair_cc_analytic %>% 
  dummy_cols(select_columns = c("sex_dv_t0",
                                "non_white_t0",
                                "marital_status_t0",
                                "hiqual_dv_t0",
                                "gor_dv_t0",
                                "sic2007_section_lab_t0",
                                "soc2000_major_group_title_t0",
                                "jbft_dv_t0",
                                "small_firm_t0",
                                "emp_contract_t0",
                                "broken_emp_t0",
                                "j2has_dv_t0",
                                "rel_pov_t0",
                                "health_t0",
                                "srh_bin_t0",
                                "ghq_case4_t0",
                                "marital_status_t1",
                                "health_t1",
                                "srh_bin_t1",
                                "ghq_case4_t1"))

pair_cc_analytic <- pair_cc_analytic %>% janitor::clean_names()

#### create age-stratified df's ------------------------------------------------
younger_df <- pair_cc_analytic %>% filter(age_dv_t0<=median(age_dv_t0))
older_df <- pair_cc_analytic %>% filter(age_dv_t0>median(age_dv_t0))


################################################################################
#####                               functions                              #####
################################################################################

#### propensity score model - mlm ----------------------------------------------
ps_model_mlm <- function(data = pair_cc_ps, outcome){
  glmmTMB::glmmTMB(outcome ~
                     sex_dv_t0 +
                     # age_dv_t0 +
                     non_white_t0 +
                     #        marital_status_t0_married_civil_partnership +
                     marital_status_t0_divorced_separated_widowed +
                     marital_status_t0_single +
                     #        hiqual_dv_t0_degree +
                     hiqual_dv_t0_other_higher_degree +
                     hiqual_dv_t0_a_level_etc +
                     hiqual_dv_t0_gcse_etc +
                     hiqual_dv_t0_other_qualification +
                     hiqual_dv_t0_no_qualification +
                     #  gor_dv_t0_east_midlands +
                     gor_dv_t0_east_of_england +
                     gor_dv_t0_london +
                     gor_dv_t0_north_east +
                     gor_dv_t0_north_west +
                     gor_dv_t0_northern_ireland +
                     gor_dv_t0_scotland +
                     gor_dv_t0_south_east +
                     gor_dv_t0_south_west +
                     gor_dv_t0_wales +
                     gor_dv_t0_west_midlands +
                     gor_dv_t0_yorkshire_and_the_humber +
                     #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                     sic2007_section_lab_t0_administrative_and_support_service_activities +
                     sic2007_section_lab_t0_construction +
                     sic2007_section_lab_t0_education +
                     sic2007_section_lab_t0_human_health_and_social_work_activities +
                     sic2007_section_lab_t0_manufacturing +
                     sic2007_section_lab_t0_other_industry +
                     sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                     sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                     sic2007_section_lab_t0_transportation_and_storage +
                     sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                     #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                     soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                     soc2000_major_group_title_t0_elementary_occupations +
                     soc2000_major_group_title_t0_managers_and_senior_officials +
                     soc2000_major_group_title_t0_personal_service_occupations +
                     soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                     soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                     soc2000_major_group_title_t0_science_and_technology_professionals +
                     soc2000_major_group_title_t0_skilled_trades_occupations +
                     jbft_dv_t0 +
                     small_firm_t0 +
                     emp_contract_t0 +
                     #  broken_emp_t0_broken_employment +
                     broken_emp_t0_no_employment_spells +
                     broken_emp_t0_unbroken_employment +
                     j2has_dv_t0 +
                     rel_pov_t0 +
                     health_t0 +
                     srh_bin_t0 +
                     ghq_case4_t0 +
                     sf12mcs_dv_t0 +
                     sf12pcs_dv_t0 +                
                     (1|pidp),
                   family = binomial(link="logit"),
                   data = data)#, control=glmerControl(optimizer="bobyqa",
  #                     optCtrl=list(maxfun=2e5)))
  
  
}

################################################################################
#####               Unweighted study characteristics tables                #####
################################################################################

younger_df <- droplevels(younger_df)
older_df <- droplevels(older_df)


#### younger --------------------------------------------------------------------
table_one_young <- tableone::CreateTableOne(vars = cov_vector3, strata = "exposure1",
                                        data = younger_df,
                                        test = FALSE)

table_one_young_sav <- print(table_one_young, showAllLevels = TRUE, smd = TRUE,
                         nonnormal = nonnorm_vec,
                         formatOptions = list(big.mark = ","))

# Count covariates with important imbalance
table_one_young_smd <- data.frame(ExtractSmd(table_one_young))
table_one_young_smd <- table_one_young_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_young_smd, "./working_data/cc/subgroup/age/table_one_young_unmatched_smd.csv")

#### older --------------------------------------------------------------------
table_one_old <- tableone::CreateTableOne(vars = cov_vector3, strata = "exposure1",
                                        data = older_df,
                                        test = FALSE)

table_one_old_sav <- print(table_one_old, showAllLevels = TRUE, smd = TRUE,
                         nonnormal = nonnorm_vec,
                         formatOptions = list(big.mark = ","))

# Count covariates with important imbalance
table_one_old_smd <- data.frame(ExtractSmd(table_one_old))
table_one_old_smd <- table_one_old_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_old_smd, "./working_data/cc/subgroup/age/table_one_old_unmatched_smd.csv")

################################################################################
#####                     Unemployment at t1 PS model                      #####
################################################################################

#### younger --------------------------------------------------------------------
### call the function (and benchmark time - takes a while to run for MLM)
start_time <- Sys.time()
ps_mod_young_MLM <- ps_model_mlm(data = younger_df, outcome = younger_df$exposure1)
end_time <- Sys.time()
end_time - start_time

### summary of model
summary(ps_mod_young_MLM)


### predicted probability of being assigned to exposed group
younger_df$ps_exp1 <- predict(ps_mod_young_MLM, type = "response")
summary(younger_df$ps_exp1)

### predicted probability of being assigned to unexposed group
younger_df$ps_noexp1 <- 1-younger_df$ps_exp1
summary(younger_df$ps_noexp1)

### predicted probability of being assigned to actual exposure status
younger_df$ps_assign <- NA
younger_df$ps_assign[younger_df$exposure1=="exposed (employed at t1)"] <- younger_df$ps_exp1[younger_df$exposure1=="exposed (employed at t1)"]
younger_df$ps_assign[younger_df$exposure1=="unexposed"] <- younger_df$ps_noexp1[younger_df$exposure1=="unexposed"]

### smaller of ps_exp1 and ps_noexp1 for matching weight
younger_df$ps_min <- pmin(younger_df$ps_exp1, younger_df$ps_noexp1)
summary(younger_df$ps_min)


#### older ---------------------------------------------------------------------
### call the function (and benchmark time - takes a while to run for MLM)
start_time <- Sys.time()
ps_mod_old_MLM <- ps_model_mlm(data = older_df, outcome = older_df$exposure1)
end_time <- Sys.time()
end_time - start_time

### summary of model
summary(ps_mod_old_MLM)


### predicted probability of being assigned to exposed group
older_df$ps_exp1 <- predict(ps_mod_old_MLM, type = "response")
summary(older_df$ps_exp1)

### predicted probability of being assigned to unexposed group
older_df$ps_noexp1 <- 1-older_df$ps_exp1
summary(older_df$ps_noexp1)

### predicted probability of being assigned to actual exposure status
older_df$ps_assign <- NA
older_df$ps_assign[older_df$exposure1=="exposed (employed at t1)"] <- older_df$ps_exp1[older_df$exposure1=="exposed (employed at t1)"]
older_df$ps_assign[older_df$exposure1=="unexposed"] <- older_df$ps_noexp1[older_df$exposure1=="unexposed"]

### smaller of ps_exp1 and ps_noexp1 for matching weight
older_df$ps_min <- pmin(older_df$ps_exp1, older_df$ps_noexp1)
summary(older_df$ps_min)


################################################################################
#####                       propensity score matching                      #####
################################################################################

#### younger -------------------------------------------------------------------
### match propensity scores using Matching package ----
list_match_young <- Match(Tr = (younger_df$exposure1=="exposed (employed at t1)"),
                      # logit of PS/1-PS
                      X = log(younger_df$ps_exp1/younger_df$ps_noexp1),
                      ## 1:1 matching ratio
                      M = 1,
                      ## caliper = 0.2 * SD(logit(PS))
                      caliper  = 0.2,
                      replace  = FALSE,
                      ties     = TRUE,
                      version  = "fast")

### extract matched data -----------------
young_matched <- younger_df[unlist(list_match_young[c("index.treated","index.control")]), ]

### matched Table One with SMD ---------
table_one_young_matched <- CreateTableOne(vars = cov_vector3,
                                      strata = "exposure1",
                                      data = young_matched,
                                      test = FALSE)

table_one_young_matched_sav <- print(table_one_young_matched,  smd = TRUE)

write.csv(table_one_young_matched_sav, "./output/matched_descriptives/subgroups/age/table_one_young_matched.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(table_one_young_matched) > 0.1))

table_one_young_matched_smd <- data.frame(ExtractSmd(table_one_young_matched))
table_one_young_matched_smd <- table_one_young_matched_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_young_matched_smd, "./working_data/cc/subgroup/age/table_one_young_matched_smd.csv")

#### older ---------------------------------------------------------------------
### match propensity scores using Matching package ----
list_match_old <- Match(Tr = (older_df$exposure1=="exposed (employed at t1)"),
                      # logit of PS/1-PS
                      X = log(older_df$ps_exp1/older_df$ps_noexp1),
                      ## 1:1 matching ratio
                      M = 1,
                      ## caliper = 0.2 * SD(logit(PS))
                      caliper  = 0.2,
                      replace  = FALSE,
                      ties     = TRUE,
                      version  = "fast")

### extract matched data -----------------
old_matched <- older_df[unlist(list_match_old[c("index.treated","index.control")]), ]

### matched Table One with SMD ---------
table_one_old_matched <- CreateTableOne(vars = cov_vector3,
                                      strata = "exposure1",
                                      data = old_matched,
                                      test = FALSE)

table_one_old_matched_sav <- print(table_one_old_matched,  smd = TRUE)

write.csv(table_one_old_matched_sav, "./output/matched_descriptives/subgroups/age/table_one_old_matched.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(table_one_old_matched) > 0.1))

table_one_old_matched_smd <- data.frame(ExtractSmd(table_one_old_matched))
table_one_old_matched_smd <- table_one_old_matched_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_old_matched_smd, "./working_data/cc/subgroup/age/table_one_old_matched_smd.csv")

################################################################################
#####                      IPTW using MatchIt package                      #####
################################################################################

#### younger --------------------------------------------------------------------
young_iptw <- younger_df %>% 
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           1,0))

young_iptw$srh_bin_t1 <- as.character(young_iptw$srh_bin_t1)
young_iptw$ghq_case4_t1 <- as.character(young_iptw$ghq_case4_t1)


start_time <- Sys.time()
young_iptw_mod <- matchit(exp1_bin ~
                        sex_dv_t0 +
                        non_white_t0 +
                        #        marital_status_t0_married_civil_partnership +
                        marital_status_t0_divorced_separated_widowed +
                        marital_status_t0_single +
                        #        hiqual_dv_t0_degree +
                        hiqual_dv_t0_other_higher_degree +
                        hiqual_dv_t0_a_level_etc +
                        hiqual_dv_t0_gcse_etc +
                        hiqual_dv_t0_other_qualification +
                        hiqual_dv_t0_no_qualification +
                        #  gor_dv_t0_east_midlands +
                        gor_dv_t0_east_of_england +
                        gor_dv_t0_london +
                        gor_dv_t0_north_east +
                        gor_dv_t0_north_west +
                        gor_dv_t0_northern_ireland +
                        gor_dv_t0_scotland +
                        gor_dv_t0_south_east +
                        gor_dv_t0_south_west +
                        gor_dv_t0_wales +
                        gor_dv_t0_west_midlands +
                        gor_dv_t0_yorkshire_and_the_humber +
                        #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                        sic2007_section_lab_t0_administrative_and_support_service_activities +
                        sic2007_section_lab_t0_construction +
                        sic2007_section_lab_t0_education +
                        sic2007_section_lab_t0_human_health_and_social_work_activities +
                        sic2007_section_lab_t0_manufacturing +
                        sic2007_section_lab_t0_other_industry +
                        sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                        sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                        sic2007_section_lab_t0_transportation_and_storage +
                        sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                        #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                        soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                        soc2000_major_group_title_t0_elementary_occupations +
                        soc2000_major_group_title_t0_managers_and_senior_officials +
                        soc2000_major_group_title_t0_personal_service_occupations +
                        soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                        soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                        soc2000_major_group_title_t0_science_and_technology_professionals +
                        soc2000_major_group_title_t0_skilled_trades_occupations +
                        jbft_dv_t0 +
                        small_firm_t0 +
                        emp_contract_t0 +
                        #  broken_emp_t0_broken_employment +
                        broken_emp_t0_no_employment_spells +
                        broken_emp_t0_unbroken_employment +
                        j2has_dv_t0 +
                        rel_pov_t0 +
                        health_t0 +
                        srh_bin_t0 +
                        ghq_case4_t0 +
                        sf12mcs_dv_t0 +
                        sf12pcs_dv_t0, 
                      data = young_iptw,
                      method = "quick", # Generalized Full Matching
                      distance = "glm",
                      estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

#test_summary <- summary(test)
#write.csv(test_summary, "./output/temp_output/test_sumary.csv")

young_iptw$weights_ps <- young_iptw_mod$weights
#matchit_df$weights <- unname(matchit_df$weights)

### create weighted data
young_iptw_svy <- svydesign(ids = ~1,
                        data = young_iptw,
                        weights = ~weights_ps)


### weighted table one
table_one_young_iptw <- svyCreateTableOne(vars = cov_vector3,
                                      strata = "exposure1",
                                      data = young_iptw_svy,
                                      test = FALSE)

table_one_young_iptw_sav <- print(table_one_young_iptw, showAllLevels = TRUE, smd = TRUE,
                              nonnormal = nonnorm_vec,
                              formatOptions = list(big.mark = ","))

write.csv(table_one_young_iptw_sav, "./output/weighted_descriptives/subgroups/age/table_one_young_iptw_sav.csv")


### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(table_one_young_iptw) > 0.1))
addmargins(table(ExtractSmd(table_one_young_iptw) > 0.2))

table_one_young_iptw_smd <- data.frame(ExtractSmd(table_one_young_iptw))
table_one_young_iptw_smd <- table_one_young_iptw_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_young_iptw_smd, "./working_data/cc/subgroup/age/table_one_young_iptw_smd.csv")

#### older ----------------------------------------------------------------------
old_iptw <- older_df %>% 
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           1,0))

old_iptw$srh_bin_t1 <- as.character(old_iptw$srh_bin_t1)
old_iptw$ghq_case4_t1 <- as.character(old_iptw$ghq_case4_t1)


start_time <- Sys.time()
old_iptw_mod <- matchit(exp1_bin ~
                        sex_dv_t0 +
                        non_white_t0 +
                        #        marital_status_t0_married_civil_partnership +
                        marital_status_t0_divorced_separated_widowed +
                        marital_status_t0_single +
                        #        hiqual_dv_t0_degree +
                        hiqual_dv_t0_other_higher_degree +
                        hiqual_dv_t0_a_level_etc +
                        hiqual_dv_t0_gcse_etc +
                        hiqual_dv_t0_other_qualification +
                        hiqual_dv_t0_no_qualification +
                        #  gor_dv_t0_east_midlands +
                        gor_dv_t0_east_of_england +
                        gor_dv_t0_london +
                        gor_dv_t0_north_east +
                        gor_dv_t0_north_west +
                        gor_dv_t0_northern_ireland +
                        gor_dv_t0_scotland +
                        gor_dv_t0_south_east +
                        gor_dv_t0_south_west +
                        gor_dv_t0_wales +
                        gor_dv_t0_west_midlands +
                        gor_dv_t0_yorkshire_and_the_humber +
                        #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                        sic2007_section_lab_t0_administrative_and_support_service_activities +
                        sic2007_section_lab_t0_construction +
                        sic2007_section_lab_t0_education +
                        sic2007_section_lab_t0_human_health_and_social_work_activities +
                        sic2007_section_lab_t0_manufacturing +
                        sic2007_section_lab_t0_other_industry +
                        sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                        sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                        sic2007_section_lab_t0_transportation_and_storage +
                        sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                        #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                        soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                        soc2000_major_group_title_t0_elementary_occupations +
                        soc2000_major_group_title_t0_managers_and_senior_officials +
                        soc2000_major_group_title_t0_personal_service_occupations +
                        soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                        soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                        soc2000_major_group_title_t0_science_and_technology_professionals +
                        soc2000_major_group_title_t0_skilled_trades_occupations +
                        jbft_dv_t0 +
                        small_firm_t0 +
                        emp_contract_t0 +
                        #  broken_emp_t0_broken_employment +
                        broken_emp_t0_no_employment_spells +
                        broken_emp_t0_unbroken_employment +
                        j2has_dv_t0 +
                        rel_pov_t0 +
                        health_t0 +
                        srh_bin_t0 +
                        ghq_case4_t0 +
                        sf12mcs_dv_t0 +
                        sf12pcs_dv_t0, 
                      data = old_iptw,
                      method = "quick", # Generalized Full Matching
                      distance = "glm",
                      estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

#test_summary <- summary(test)
#write.csv(test_summary, "./output/temp_output/test_sumary.csv")

old_iptw$weights_ps <- old_iptw_mod$weights
#matchit_df$weights <- unname(matchit_df$weights)

### create weighted data
old_iptw_svy <- svydesign(ids = ~1,
                        data = old_iptw,
                        weights = ~weights_ps)


### weighted table one
table_one_old_iptw <- svyCreateTableOne(vars = cov_vector3,
                                      strata = "exposure1",
                                      data = old_iptw_svy,
                                      test = FALSE)

table_one_old_iptw_sav <- print(table_one_old_iptw, showAllLevels = TRUE, smd = TRUE,
                              #                              nonnormal = nonnorm_vec,
                              formatOptions = list(big.mark = ","))

write.csv(table_one_old_iptw_sav, "./output/weighted_descriptives/subgroups/age/table_one_old_iptw_sav.csv")


### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(table_one_old_iptw) > 0.1))
addmargins(table(ExtractSmd(table_one_old_iptw) > 0.2))

table_one_old_iptw_smd <- data.frame(ExtractSmd(table_one_old_iptw))
table_one_old_iptw_smd <- table_one_old_iptw_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_old_iptw_smd, "./working_data/cc/subgroup/age/table_one_old_iptw_smd.csv")

################################################################################
#####                               outcomes                              ######
################################################################################

# double robust MSM

#### younger --------------------------------------------------------------------

### SF-12 PCS -----------------------
start_time <- Sys.time()
dr_iptw_pcs_young_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         sf12pcs_dv_t0 +
                                # interaction terms
                                #sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                # age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              weights = weights_ps,
                              data = young_iptw)
end_time <- Sys.time()
end_time-start_time

dr_iptw_pcs_young <- summary(dr_iptw_pcs_young_mod)

#fixef(dr_iptw_pcs_glmmTMB_mod)


dr_iptw_pcs_young_df <- tidy(dr_iptw_pcs_young_mod)

## confidence intervals
dr_iptw_pcs_young_df_ci <- data.frame(confint(dr_iptw_pcs_young_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_pcs_young_df_ci <- cbind(rownames(dr_iptw_pcs_young_df_ci),dr_iptw_pcs_young_df_ci, row.names=NULL)

dr_iptw_pcs_young_df_ci <- dr_iptw_pcs_young_df_ci %>% 
  rename(term = `rownames(dr_iptw_pcs_young_df_ci)`)

## join dfs together
dr_iptw_pcs_young_df <- dr_iptw_pcs_young_df %>% 
  left_join(dr_iptw_pcs_young_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------
start_time <- Sys.time()
dr_iptw_mcs_young_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                # age_dv_t0 +
                                # age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         sf12mcs_dv_t0 +
                                # interaction terms
                                #sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                # age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              weights = weights_ps,
                              data = young_iptw)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_mcs_young_mod)

diagnose(dr_iptw_mcs_young_mod)

dr_iptw_mcs_young_df <- tidy(dr_iptw_mcs_young_mod)

## confidence intervals
dr_iptw_mcs_young_df_ci <- data.frame(confint(dr_iptw_mcs_young_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_mcs_young_df_ci <- cbind(rownames(dr_iptw_mcs_young_df_ci),dr_iptw_mcs_young_df_ci, row.names=NULL)

dr_iptw_mcs_young_df_ci <- dr_iptw_mcs_young_df_ci %>% 
  rename(term = `rownames(dr_iptw_mcs_young_df_ci)`)

## join dfs together
dr_iptw_mcs_young_df <- dr_iptw_mcs_young_df %>% 
  left_join(dr_iptw_mcs_young_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

young_iptw$srh_bin_t1 <- factor(young_iptw$srh_bin_t1,
                            levels = c("excellent/very good", 
                                       "good/fair/poor"))


start_time <- Sys.time()
dr_iptw_srh_young_mod <- glmmTMB( srh_bin_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         srh_bin_t0 +
                                # interaction terms
                                #sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                # age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              family = binomial(link="logit"),
                              weights = weights_ps,
                              data = young_iptw)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_srh_young_mod)
diagnose(dr_iptw_srh_young_mod)

dr_iptw_srh_young_df <- tidy(dr_iptw_srh_young_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_srh_young_df_ci <- data.frame(confint(dr_iptw_srh_young_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_srh_young_df_ci <- cbind(rownames(dr_iptw_srh_young_df_ci),dr_iptw_srh_young_df_ci, row.names=NULL)

dr_iptw_srh_young_df_ci <- dr_iptw_srh_young_df_ci %>% 
  rename(term = `rownames(dr_iptw_srh_young_df_ci)`)

## join dfs together
dr_iptw_srh_young_df <- dr_iptw_srh_young_df %>% 
  left_join(dr_iptw_srh_young_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------

young_iptw$ghq_case4_t1 <- factor(young_iptw$ghq_case4_t1,
                              levels = c("0-3", "4 or more"))


start_time <- Sys.time()
dr_iptw_ghq_young_mod <- glmmTMB( ghq_case4_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         ghq_case4_t0 +
                                # interaction terms
                                #                                      sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                #age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              family = binomial(link="logit"),
                              weights = weights_ps,
                              data = young_iptw)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_ghq_young_mod)
diagnose(dr_iptw_ghq_young_mod)

dr_iptw_ghq_young_df <- tidy(dr_iptw_ghq_young_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_ghq_young_df_ci <- data.frame(confint(dr_iptw_ghq_young_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_ghq_young_df_ci <- cbind(rownames(dr_iptw_ghq_young_df_ci),dr_iptw_ghq_young_df_ci, row.names=NULL)

dr_iptw_ghq_young_df_ci <- dr_iptw_ghq_young_df_ci %>% 
  rename(term = `rownames(dr_iptw_ghq_young_df_ci)`)

## join dfs together
dr_iptw_ghq_young_df <- dr_iptw_ghq_young_df %>% 
  left_join(dr_iptw_ghq_young_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### create single summary df for iptw double-robust cc outcomes ------------

dr_iptw_young_df <- dr_iptw_pcs_young_df %>% 
  bind_rows(dr_iptw_mcs_young_df,
            dr_iptw_srh_young_df,
            dr_iptw_ghq_young_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, effect, est_type, estimate, std.error, p.value, lci, uci)

write.csv(dr_iptw_young_df, "./output/weighted_outcomes/cc/sub_groups/age/dr_iptw_young_df.csv")


diagnose(dr_iptw_pcs_young_mod)
diagnose(dr_iptw_mcs_young_mod)
diagnose(dr_iptw_srh_young_mod)
diagnose(dr_iptw_ghq_young_mod)




#### older ----------------------------------------------------------------------

### SF-12 PCS -----------------------
start_time <- Sys.time()
dr_iptw_pcs_old_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         sf12pcs_dv_t0 +
                                # interaction terms
                                #                                      sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                #age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              weights = weights_ps,
                              data = old_iptw)
end_time <- Sys.time()
end_time-start_time

dr_iptw_pcs_old <- summary(dr_iptw_pcs_old_mod)

#fixef(dr_iptw_pcs_glmmTMB_mod)


dr_iptw_pcs_old_df <- tidy(dr_iptw_pcs_old_mod)

## confidence intervals
dr_iptw_pcs_old_df_ci <- data.frame(confint(dr_iptw_pcs_old_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_pcs_old_df_ci <- cbind(rownames(dr_iptw_pcs_old_df_ci),dr_iptw_pcs_old_df_ci, row.names=NULL)

dr_iptw_pcs_old_df_ci <- dr_iptw_pcs_old_df_ci %>% 
  rename(term = `rownames(dr_iptw_pcs_old_df_ci)`)

## join dfs together
dr_iptw_pcs_old_df <- dr_iptw_pcs_old_df %>% 
  left_join(dr_iptw_pcs_old_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------
start_time <- Sys.time()
dr_iptw_mcs_old_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         sf12mcs_dv_t0 +
                                # interaction terms
                                #                                     sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                #age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              weights = weights_ps,
                              data = old_iptw)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_mcs_old_mod)

diagnose(dr_iptw_mcs_old_mod)

dr_iptw_mcs_old_df <- tidy(dr_iptw_mcs_old_mod)

## confidence intervals
dr_iptw_mcs_old_df_ci <- data.frame(confint(dr_iptw_mcs_old_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_mcs_old_df_ci <- cbind(rownames(dr_iptw_mcs_old_df_ci),dr_iptw_mcs_old_df_ci, row.names=NULL)

dr_iptw_mcs_old_df_ci <- dr_iptw_mcs_old_df_ci %>% 
  rename(term = `rownames(dr_iptw_mcs_old_df_ci)`)

## join dfs together
dr_iptw_mcs_old_df <- dr_iptw_mcs_old_df %>% 
  left_join(dr_iptw_mcs_old_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

old_iptw$srh_bin_t1 <- factor(old_iptw$srh_bin_t1,
                            levels = c("excellent/very good", 
                                       "good/fair/poor"))


start_time <- Sys.time()
dr_iptw_srh_old_mod <- glmmTMB( srh_bin_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         srh_bin_t0 +
                                # interaction terms
                                #                                      sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                #age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              family = binomial(link="logit"),
                              weights = weights_ps,
                              data = old_iptw)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_srh_old_mod)
diagnose(dr_iptw_srh_old_mod)

dr_iptw_srh_old_df <- tidy(dr_iptw_srh_old_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_srh_old_df_ci <- data.frame(confint(dr_iptw_srh_old_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_srh_old_df_ci <- cbind(rownames(dr_iptw_srh_old_df_ci),dr_iptw_srh_old_df_ci, row.names=NULL)

dr_iptw_srh_old_df_ci <- dr_iptw_srh_old_df_ci %>% 
  rename(term = `rownames(dr_iptw_srh_old_df_ci)`)

## join dfs together
dr_iptw_srh_old_df <- dr_iptw_srh_old_df %>% 
  left_join(dr_iptw_srh_old_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------

old_iptw$ghq_case4_t1 <- factor(old_iptw$ghq_case4_t1,
                              levels = c("0-3", "4 or more"))


start_time <- Sys.time()
dr_iptw_ghq_old_mod <- glmmTMB( ghq_case4_t1 ~
                                exposure1 +
                                sex_dv_t0 +
                                #age_dv_t0 +
                                #age_dv_t1 +
                                non_white_t0 +
                                #        marital_status_t0_married_civil_partnership +
                                marital_status_t0_divorced_separated_widowed +
                                marital_status_t0_single +
                                #        marital_status_t1_married_civil_partnership +
                                marital_status_t1_divorced_separated_widowed +
                                marital_status_t1_single +
                                #        hiqual_dv_t0_degree +
                                hiqual_dv_t0_other_higher_degree +
                                hiqual_dv_t0_a_level_etc +
                                hiqual_dv_t0_gcse_etc +
                                hiqual_dv_t0_other_qualification +
                                hiqual_dv_t0_no_qualification +
                                #  gor_dv_t0_east_midlands +
                                gor_dv_t0_east_of_england +
                                gor_dv_t0_london +
                                gor_dv_t0_north_east +
                                gor_dv_t0_north_west +
                                gor_dv_t0_northern_ireland +
                                gor_dv_t0_scotland +
                                gor_dv_t0_south_east +
                                gor_dv_t0_south_west +
                                gor_dv_t0_wales +
                                gor_dv_t0_west_midlands +
                                gor_dv_t0_yorkshire_and_the_humber +
                                #                                         #  gor_dv_t1_east_midlands +
                                #                                         gor_dv_t1_east_of_england +
                                #                                         gor_dv_t1_london +
                                #                                         gor_dv_t1_north_east +
                                #                                         gor_dv_t1_north_west +
                                #                                         gor_dv_t1_northern_ireland +
                                #                                         gor_dv_t1_scotland +
                                #                                         gor_dv_t1_south_east +
                                #                                         gor_dv_t1_south_west +
                                #                                         gor_dv_t1_wales +
                                #                                         gor_dv_t1_west_midlands +
                              #                                         gor_dv_t1_yorkshire_and_the_humber +
                              #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
                              sic2007_section_lab_t0_administrative_and_support_service_activities +
                                sic2007_section_lab_t0_construction +
                                sic2007_section_lab_t0_education +
                                sic2007_section_lab_t0_human_health_and_social_work_activities +
                                sic2007_section_lab_t0_manufacturing +
                                sic2007_section_lab_t0_other_industry +
                                sic2007_section_lab_t0_professional_scientific_and_technical_activities +
                                sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
                                sic2007_section_lab_t0_transportation_and_storage +
                                sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
                                #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
                                soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
                                soc2000_major_group_title_t0_elementary_occupations +
                                soc2000_major_group_title_t0_managers_and_senior_officials +
                                soc2000_major_group_title_t0_personal_service_occupations +
                                soc2000_major_group_title_t0_process_plant_and_machine_operatives +
                                soc2000_major_group_title_t0_sales_and_customer_service_occupations +
                                soc2000_major_group_title_t0_science_and_technology_professionals +
                                soc2000_major_group_title_t0_skilled_trades_occupations +
                                jbft_dv_t0 +
                                small_firm_t0 +
                                emp_contract_t0 +
                                #  broken_emp_t0_broken_employment +
                                broken_emp_t0_no_employment_spells +
                                broken_emp_t0_unbroken_employment +
                                j2has_dv_t0 +
                                rel_pov_t0 +
                                #                                        rel_pov_t1 +
                                health_t0 +
                                health_t1 +
                                #                                         ghq_case4_t0 +
                                # interaction terms
                                #                                      sex_dv_t0*age_dv_t0 +
                                sex_dv_t0*rel_pov_t0 +
                                #age_dv_t0*rel_pov_t0 +
                                (1|pidp),
                              family = binomial(link="logit"),
                              weights = weights_ps,
                              data = old_iptw)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_ghq_old_mod)
diagnose(dr_iptw_ghq_old_mod)

dr_iptw_ghq_old_df <- tidy(dr_iptw_ghq_old_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_ghq_old_df_ci <- data.frame(confint(dr_iptw_ghq_old_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_ghq_old_df_ci <- cbind(rownames(dr_iptw_ghq_old_df_ci),dr_iptw_ghq_old_df_ci, row.names=NULL)

dr_iptw_ghq_old_df_ci <- dr_iptw_ghq_old_df_ci %>% 
  rename(term = `rownames(dr_iptw_ghq_old_df_ci)`)

## join dfs together
dr_iptw_ghq_old_df <- dr_iptw_ghq_old_df %>% 
  left_join(dr_iptw_ghq_old_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### create single summary df for iptw double-robust cc outcomes ------------

dr_iptw_old_df <- dr_iptw_pcs_old_df %>% 
  bind_rows(dr_iptw_mcs_old_df,
            dr_iptw_srh_old_df,
            dr_iptw_ghq_old_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, effect, est_type, estimate, std.error, p.value, lci, uci)

write.csv(dr_iptw_old_df, "./output/weighted_outcomes/cc/sub_groups/age/dr_iptw_old_df.csv")


diagnose(dr_iptw_pcs_old_mod)
diagnose(dr_iptw_mcs_old_mod)
diagnose(dr_iptw_srh_old_mod)
diagnose(dr_iptw_ghq_old_mod)
