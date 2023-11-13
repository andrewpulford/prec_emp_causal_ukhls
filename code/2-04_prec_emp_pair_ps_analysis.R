################################################################################

# Precarious employment and health - Understanding Society
# 3-05 - Paired propensity score matched complete case outcome analysis for risk  
# of job loss 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) complete case outcome analysis for propensity score matched data



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

#### load propensity matched df ------------------------------------------------
df_matched <- readRDS("working_data/pair_cc_matched.rds") # PS matched data

#### load analytic df
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")

pair_cc_analytic$pidp <- factor(pair_cc_analytic$pidp)
pair_cc_analytic$sf12mcs_dv_t0 <- as.double(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t0 <- as.double(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t1 <- as.double(pair_cc_analytic$sf12mcs_dv_t1)
pair_cc_analytic$sf12pcs_dv_t1 <- as.double(pair_cc_analytic$sf12pcs_dv_t1)

df_matched <- df_matched %>% 
  left_join(pair_cc_analytic)

df_matched <- df_matched %>% 
  dummy_cols(select_columns = c("marital_status_t1",
                                "gor_dv_t1",
                                "health_t1"
  ))

df_matched <- df_matched %>% janitor::clean_names()


### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
df_matched$srh_bin_t0 <- factor(df_matched$srh_bin_t0,
                                levels = c("excellent/very good", 
                                           "good/fair/poor"))
df_matched$srh_bin_t1 <- factor(df_matched$srh_bin_t1,
                                levels = c("excellent/very good", 
                                           "good/fair/poor"))

df_matched$ghq_case4_t0 <- factor(df_matched$ghq_case4_t0,
                                  levels = c("0-3", "4 or more"))
df_matched$ghq_case4_t1 <- factor(df_matched$ghq_case4_t1,
                                  levels = c("0-3", "4 or more"))


df_matched$exposure1 <- factor(df_matched$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))
df_matched$exposure2 <- factor(df_matched$exposure2,
                               levels = c("unexposed",
                                          "exposed (no job loss between t0 and t1"))
# error with exposure 2 to correct at some point

### convert SF-12 outcomes to numeric to allow svyglm to work
df_matched$sf12pcs_dv_t1 <- as.numeric(df_matched$sf12pcs_dv_t1)
df_matched$sf12mcs_dv_t1 <- as.numeric(df_matched$sf12mcs_dv_t1)



#### load IPTW df --------------------------------------------------------------
matchit_df <- readRDS("working_data/matchit_df.rds") # analytic df with IPTW from MatchIt package


### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
matchit_df$srh_bin_t0 <- factor(matchit_df$srh_bin_t0,
                                      levels = c("good/fair/poor", 
                                                 "excellent/very good"))
matchit_df$srh_bin_t1 <- factor(matchit_df$srh_bin_t1,
                                      levels = c("good/fair/poor", 
                                                 "excellent/very good"))

matchit_df$ghq_case4_t0 <- factor(matchit_df$ghq_case4_t0,
                                        levels = c("0-3", "4 or more"))
matchit_df$ghq_case4_t1 <- factor(matchit_df$ghq_case4_t1,
                                        levels = c("0-3", "4 or more"))


matchit_df$exposure1 <- factor(matchit_df$exposure1,
                                     levels = c("unexposed",
                                                "exposed (employed at t1)"))
matchit_df$exposure2 <- factor(matchit_df$exposure2,
                                     levels = c("unexposed",
                                                "exposed (no job loss between t0 and t1"))
# error with exposure 2

### convert SF-12 outcomes to numeric to allow svyglm to work
matchit_df$sf12pcs_dv_t1 <- as.numeric(matchit_df$sf12pcs_dv_t1)
matchit_df$sf12mcs_dv_t1 <- as.numeric(matchit_df$sf12mcs_dv_t1)

### create weighted data
matchit_svy <- svydesign(ids = ~1,
                         data = matchit_df,
                         weights = ~weights_ps)




################################################################################
#####                           descriptive tables                         #####
################################################################################

#### unweighted outcomes table -------------------------------------------------
table_outcomes_unweighted <- CreateTableOne(vars = outcome_vector2, 
                                            strata = "exposure1",
                                            data = matchit_df,
                                            test = TRUE)

table_outcomes_unweighted_sav <- print(table_outcomes_unweighted, 
                                       showAllLevels = TRUE, 
                                       smd = TRUE, 
                                       formatOptions = list(big.mark = ","))

#### PS matched outcomes table -------------------------------------------------
table_outcomes_matched <- CreateTableOne(vars = outcome_vector2, 
                                            strata = "exposure1",
                                            data = df_matched,
                                            test = TRUE)

table_outcomes_matched_sav <- print(table_outcomes_matched, 
                                       showAllLevels = TRUE, 
                                       smd = TRUE, 
                                       formatOptions = list(big.mark = ","))


#### weighted outcomes table ---------------------------------------------------
table_outcomes_IPTW <- svyCreateTableOne(vars = outcome_vector2,
                                             strata = "exposure1",
                                             data = matchit_svy,
                                             test = TRUE)

table_outcomes_IPTW_sav <- print(table_outcomes_IPTW, 
                                     showAllLevels = TRUE,  
                                     smd = TRUE)

write.csv(table_outcomes_IPTW_sav, "./output/weighted_descriptives/table_outcomes_IPTW_sav.csv")

################################################################################
#####                       matched regression models                     ######
################################################################################

# initial models prior to doubly robust estimation 

### SF-12 PCS -----------------------
matched_pcs_glmmTMB_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                      exposure1 +
           (1|pidp),
         data = df_matched)

summary(matched_pcs_glmmTMB_mod)

### SF-12 MCS -----------------------
matched_mcs_glmmTMB_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                      exposure1 +
                                      (1|pidp),
                                    data = df_matched)

summary(matched_mcs_glmmTMB_mod)

### poor self-rated health -------------------

matched_srh_glmmTMB_mod <- glmmTMB(srh_bin_t1 ~ exposure1  +
                                     (1|pidp),
                                   family = binomial(link="logit"),
                                   data = df_matched, 
                         na.action = na.omit)

summary(matched_srh_glmmTMB_mod)

### GHQ-12 caseness (4+) -------------------

matched_ghq_glmmTMB_mod <- glmmTMB(ghq_case4_t1 ~ exposure1  +
                                     (1|pidp),
                                   family = binomial(link="logit"),
                                   data = df_matched, 
                                   na.action = na.omit)

summary(matched_ghq_glmmTMB_mod)


################################################################################
#####               double robust matched regression models               ######
################################################################################



### SF-12 PCS -----------------------
start_time <- Sys.time()
dr_matched_pcs_glmmTMB_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                         exposure1 +
                                         sex_dv_t0 +
                                         age_dv_t0 +
                                         age_dv_t1 +
                                         non_white_t0 +
                                         #        marital_status_t0_married_civil_partnership +
                                         marital_status_t0_divorced_separated_widowed +
                                         marital_status_t0_single +
                                         #        marital_status_t1_married_civil_partnership +
                                         marital_status_t1_divorced_separated_widowed +
                                         marital_status_t1_single +
                                         dep_child_bin_t0 +
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
                                         health_t1 +
                                      # interaction terms
                                      sex_dv_t0*age_dv_t0 +
                                      sex_dv_t0*rel_pov_t0 +
                                      age_dv_t0*rel_pov_t0 +
                                      (1|pidp),
                                    data = df_matched)
end_time <- Sys.time()
end_time-start_time

summary(dr_matched_pcs_glmmTMB_mod)

dr_matched_pcs_df <- tidy(dr_matched_pcs_glmmTMB_mod)

## confidence intervals
dr_matched_pcs_df_ci <- data.frame(confint(dr_matched_pcs_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_matched_pcs_df_ci <- cbind(rownames(dr_matched_pcs_df_ci),dr_matched_pcs_df_ci, row.names=NULL)

dr_matched_pcs_df_ci <- dr_matched_pcs_df_ci %>% 
  rename(term = `rownames(dr_matched_pcs_df_ci)`)

## join dfs together
dr_matched_pcs_df <- dr_matched_pcs_df %>% 
  left_join(dr_matched_pcs_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------
start_time <- Sys.time()
dr_matched_mcs_glmmTMB_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                         exposure1 +
                                         sex_dv_t0 +
                                         age_dv_t0 +
                                         age_dv_t1 +
                                         non_white_t0 +
                                         #        marital_status_t0_married_civil_partnership +
                                         marital_status_t0_divorced_separated_widowed +
                                         marital_status_t0_single +
                                         #        marital_status_t1_married_civil_partnership +
                                         marital_status_t1_divorced_separated_widowed +
                                         marital_status_t1_single +
                                         dep_child_bin_t0 +
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
                                         health_t1 +
                                         sf12mcs_dv_t0 +
  # interaction terms
  sex_dv_t0*age_dv_t0 +
  sex_dv_t0*rel_pov_t0 +
  age_dv_t0*rel_pov_t0 +
                                           (1|pidp),
                                       data = df_matched)
end_time <- Sys.time()
end_time-start_time

summary(dr_matched_mcs_glmmTMB_mod)

diagnose(dr_matched_mcs_glmmTMB_mod)

dr_matched_mcs_df <- tidy(dr_matched_mcs_glmmTMB_mod)

## confidence intervals
dr_matched_mcs_df_ci <- data.frame(confint(dr_matched_mcs_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_matched_mcs_df_ci <- cbind(rownames(dr_matched_mcs_df_ci),dr_matched_mcs_df_ci, row.names=NULL)

dr_matched_mcs_df_ci <- dr_matched_mcs_df_ci %>% 
  rename(term = `rownames(dr_matched_mcs_df_ci)`)

## join dfs together
dr_matched_mcs_df <- dr_matched_mcs_df %>% 
  left_join(dr_matched_mcs_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

start_time <- Sys.time()
dr_matched_srh_glmmTMB_mod <- glmmTMB( srh_bin_t1 ~
                                         exposure1 +
                                         sex_dv_t0 +
                                         age_dv_t0 +
                                         age_dv_t1 +
                                         non_white_t0 +
                                         #        marital_status_t0_married_civil_partnership +
                                         marital_status_t0_divorced_separated_widowed +
                                         marital_status_t0_single +
                                         #        marital_status_t1_married_civil_partnership +
                                         marital_status_t1_divorced_separated_widowed +
                                         marital_status_t1_single +
                                         dep_child_bin_t0 +
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
                                         health_t1 +
    # interaction terms
  sex_dv_t0*age_dv_t0 +
  sex_dv_t0*rel_pov_t0 +
  age_dv_t0*rel_pov_t0 +
  (1|pidp),
                                       family = binomial(link="logit"),
                                       data = df_matched)
end_time <- Sys.time()
end_time-start_time

summary(dr_matched_srh_glmmTMB_mod)

dr_matched_srh_df <- tidy(dr_matched_srh_glmmTMB_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_matched_srh_df_ci <- data.frame(confint(dr_matched_srh_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_matched_srh_df_ci <- cbind(rownames(dr_matched_srh_df_ci),dr_matched_srh_df_ci, row.names=NULL)

dr_matched_srh_df_ci <- dr_matched_srh_df_ci %>% 
  rename(term = `rownames(dr_matched_srh_df_ci)`)

## join dfs together
dr_matched_srh_df <- dr_matched_srh_df %>% 
  left_join(dr_matched_srh_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------
start_time <- Sys.time()
dr_matched_ghq_glmmTMB_mod <- glmmTMB( ghq_case4_t1 ~
                                         exposure1 +
                                         sex_dv_t0 +
                                         age_dv_t0 +
                                         age_dv_t1 +
                                         non_white_t0 +
                                         #        marital_status_t0_married_civil_partnership +
                                         marital_status_t0_divorced_separated_widowed +
                                         marital_status_t0_single +
                                         #        marital_status_t1_married_civil_partnership +
                                         marital_status_t1_divorced_separated_widowed +
                                         marital_status_t1_single +
                                         dep_child_bin_t0 +
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
                                         health_t1 +
                                         ghq_case4_t0 +
   # interaction terms
   sex_dv_t0*age_dv_t0 +
   sex_dv_t0*rel_pov_t0 +
   age_dv_t0*rel_pov_t0 +
   (1|pidp),
                                       family = binomial(link="logit"),
                                       data = df_matched)
end_time <- Sys.time()
end_time-start_time

summary(dr_matched_ghq_glmmTMB_mod)

dr_matched_ghq_df <- tidy(dr_matched_ghq_glmmTMB_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_matched_ghq_df_ci <- data.frame(confint(dr_matched_ghq_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_matched_ghq_df_ci <- cbind(rownames(dr_matched_ghq_df_ci),dr_matched_ghq_df_ci, row.names=NULL)

dr_matched_ghq_df_ci <- dr_matched_ghq_df_ci %>% 
  rename(term = `rownames(dr_matched_ghq_df_ci)`)

## join dfs together
dr_matched_ghq_df <- dr_matched_ghq_df %>% 
  left_join(dr_matched_ghq_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### create single summary df for matched double-robust cc outcomes ------------

dr_matched_df <- dr_matched_pcs_df %>% 
  bind_rows(dr_matched_mcs_df,
            dr_matched_srh_df,
            dr_matched_ghq_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, effect, est_type, estimate, std.error, p.value, lci, uci)

write.csv(dr_matched_df, "./output/matched_outcomes/cc/dr_matched_df.csv")

diagnose(dr_matched_pcs_glmmTMB_mod)
diagnose(dr_matched_mcs_glmmTMB_mod)
diagnose(dr_matched_srh_glmmTMB_mod)
diagnose(dr_matched_ghq_glmmTMB_mod)

