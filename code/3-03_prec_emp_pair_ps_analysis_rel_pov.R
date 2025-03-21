################################################################################

# Precarious employment and health - Understanding Society
# 3-03	Complete case outcome analysis (by relative poverty) – propensity score matched data 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) complete case descriptive and outcome analysis by relative poverty for PS matched data
# (b) complete case descriptive and outcome analysis by relative poverty for IPTW data


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
library(MatchIt) # PS matching
library(WeightIt) # IPTW
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
pair_cc_analytic <- readRDS("./working_data/cc/pair_cc_analytic.rds")


### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
pair_cc_analytic$srh_bin_t0 <- factor(pair_cc_analytic$srh_bin_t0,
                                      levels = c("excellent/very good/good", 
                                                 "fair/poor"))
pair_cc_analytic$srh_bin_t1 <- factor(pair_cc_analytic$srh_bin_t1,
                                      levels = c("excellent/very good/good", 
                                                 "fair/poor"))

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
                                "gor_dv_t0",
                                "sic2007_section_lab_t0",
                                "soc2000_major_group_title_t0",
                                "jbft_dv_t0",
                                "small_firm_t0",
                                "emp_contract_t0",
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
rel_pov_df <- pair_cc_analytic %>% filter(rel_pov_t0=="relative poverty")
not_pov_df <- pair_cc_analytic %>% filter(rel_pov_t0=="not in relative poverty")


################################################################################
#####                               functions                              #####
################################################################################

#### propensity score model - mlm ----------------------------------------------
ps_model_mlm <- function(data = pair_cc_ps, outcome){
  glmmTMB::glmmTMB(outcome ~
                     sex_dv_t0 +
                     age_dv_t0 +
                     non_white_t0 +
                     #        marital_status_t0_married_civil_partnership +
                     marital_status_t0_divorced_separated_widowed +
                     marital_status_t0_single +
                     dep_child_bin_t0 +
                     degree_bin_t0 +
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
                     broken_emp_t0 +
                     j2has_dv_t0 +
                     #rel_pov_t0 +
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

rel_pov_df <- droplevels(rel_pov_df)
not_pov_df <- droplevels(not_pov_df)


#### rel_pov --------------------------------------------------------------------
table_one_rel_pov <- tableone::CreateTableOne(vars = cov_vector3, strata = "exposure1",
                                            data = rel_pov_df,
                                            test = FALSE)

table_one_rel_pov_sav <- print(table_one_rel_pov, showAllLevels = TRUE, smd = TRUE,
                             nonnormal = nonnorm_vec,
                             formatOptions = list(big.mark = ","))

# Count covariates with important imbalance
table_one_rel_pov_smd <- data.frame(ExtractSmd(table_one_rel_pov))
table_one_rel_pov_smd <- table_one_rel_pov_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_rel_pov_smd, "./working_data/cc/subgroup/rel_poverty/table_one_rel_pov_unmatched_smd.csv")

#### not_pov --------------------------------------------------------------------
table_one_not_pov <- tableone::CreateTableOne(vars = cov_vector3, strata = "exposure1",
                                          data = not_pov_df,
                                          test = FALSE)

table_one_not_pov_sav <- print(table_one_not_pov, showAllLevels = TRUE, smd = TRUE,
                           nonnormal = nonnorm_vec,
                           formatOptions = list(big.mark = ","))

# Count covariates with important imbalance
table_one_not_pov_smd <- data.frame(ExtractSmd(table_one_not_pov))
table_one_not_pov_smd <- table_one_not_pov_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_not_pov_smd, "./working_data/cc/subgroup/rel_poverty/table_one_not_pov_unmatched_smd.csv")

################################################################################
#####                     Unemployment at t1 PS model                      #####
################################################################################

#### rel_pov --------------------------------------------------------------------
### call the function (and benchmark time - takes a while to run for MLM)
start_time <- Sys.time()
ps_mod_rel_pov_MLM <- ps_model_mlm(data = rel_pov_df, outcome = rel_pov_df$exposure1)
end_time <- Sys.time()
end_time - start_time

### summary of model
summary(ps_mod_rel_pov_MLM)


#### not_pov ---------------------------------------------------------------------
### call the function (and benchmark time - takes a while to run for MLM)
start_time <- Sys.time()
ps_mod_not_pov_MLM <- ps_model_mlm(data = not_pov_df, outcome = not_pov_df$exposure1)
end_time <- Sys.time()
end_time - start_time

### summary of model
summary(ps_mod_not_pov_MLM)



################################################################################
#####                       propensity score matching                      #####
################################################################################

#### rel_pov -------------------------------------------------------------------
rel_pov_matchit_df <- rel_pov_df  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed to allow 3:1 ratio matching

### convert SF-12 outcomes to numeric to allow svyglm to work
rel_pov_matchit_df$sf12pcs_dv_t0 <- as.numeric(rel_pov_matchit_df$sf12pcs_dv_t0)
rel_pov_matchit_df$sf12mcs_dv_t0 <- as.numeric(rel_pov_matchit_df$sf12mcs_dv_t0)
rel_pov_matchit_df$sf12pcs_dv_t1 <- as.numeric(rel_pov_matchit_df$sf12pcs_dv_t1)
rel_pov_matchit_df$sf12mcs_dv_t1 <- as.numeric(rel_pov_matchit_df$sf12mcs_dv_t1)

rel_pov_matchit_df$srh_bin_t1 <- as.character(rel_pov_matchit_df$srh_bin_t1)
rel_pov_matchit_df$ghq_case4_t1 <- as.character(rel_pov_matchit_df$ghq_case4_t1)


start_time <- Sys.time()
rel_pov_matchit_mod <- matchit(exp1_bin ~
                                 sex_dv_t0 +
                                 age_dv_t0 +
                                 non_white_t0 +
                                 #        marital_status_t0_married_civil_partnership +
                                 marital_status_t0_divorced_separated_widowed +
                                 marital_status_t0_single +
                                 dep_child_bin_t0 +
                                 degree_bin_t0 +
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
                                 broken_emp_t0 +
                                 j2has_dv_t0 +
                                 #rel_pov_t0 +
                                 health_t0 +
                                 srh_bin_t0 +
                                 ghq_case4_t0 +
                                 sf12mcs_dv_t0 +
                                 sf12pcs_dv_t0, 
                               data = rel_pov_matchit_df,
                               method = "nearest", 
                               ratio=3,
                               distance = "glm",
                               estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

rel_pov_matchit_df$weights_ps <- rel_pov_matchit_mod$weights

summary(rel_pov_matchit_df$weights_ps)

rel_pov_matchit_df <- rel_pov_matchit_df %>% filter(weights_ps!=0)

### create weighted data
rel_pov_matchit_df_svy <- svydesign(ids = ~1,
                                    data = rel_pov_matchit_df,
                                    weights = ~weights_ps)


### PS matched table one
rel_pov_table_one_matchit <- svyCreateTableOne(vars = cov_vector3,
                                               strata = "exposure1",
                                               data = rel_pov_matchit_df_svy,
                                               factorVars = c(catVars_vec),
                                               test = FALSE)

rel_pov_table_one_matchit_sav <- print(rel_pov_table_one_matchit, showAllLevels = TRUE, smd = TRUE,
                                       nonnormal = nonnorm_vec,
                                       factorVars = c(catVars_vec),
                                       formatOptions = list(big.mark = ","))


write.csv(rel_pov_table_one_matchit_sav, "./output/cc/matched_descriptives/subgroups/rel_poverty/rel_pov_table_one_matchit_sav.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(rel_pov_table_one_matchit) > 0.1))

rel_pov_table_one_matchit_smd <- data.frame(ExtractSmd(rel_pov_table_one_matchit))
rel_pov_table_one_matchit_smd <- rel_pov_table_one_matchit_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(rel_pov_table_one_matchit_smd, "./working_data/cc/subgroup/rel_poverty/rel_pov_table_one_matchit_smd.csv")



#### not_pov -------------------------------------------------------------------
not_pov_matchit_df <- not_pov_df  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed to allow 3:1 ratio matching

### convert SF-12 outcomes to numeric to allow svyglm to work
not_pov_matchit_df$sf12pcs_dv_t0 <- as.numeric(not_pov_matchit_df$sf12pcs_dv_t0)
not_pov_matchit_df$sf12mcs_dv_t0 <- as.numeric(not_pov_matchit_df$sf12mcs_dv_t0)
not_pov_matchit_df$sf12pcs_dv_t1 <- as.numeric(not_pov_matchit_df$sf12pcs_dv_t1)
not_pov_matchit_df$sf12mcs_dv_t1 <- as.numeric(not_pov_matchit_df$sf12mcs_dv_t1)

not_pov_matchit_df$srh_bin_t1 <- as.character(not_pov_matchit_df$srh_bin_t1)
not_pov_matchit_df$ghq_case4_t1 <- as.character(not_pov_matchit_df$ghq_case4_t1)


start_time <- Sys.time()
not_pov_matchit_mod <- matchit(exp1_bin ~
                                 sex_dv_t0 +
                                 age_dv_t0 +
                                 non_white_t0 +
                                 #        marital_status_t0_married_civil_partnership +
                                 marital_status_t0_divorced_separated_widowed +
                                 marital_status_t0_single +
                                 dep_child_bin_t0 +
                                 degree_bin_t0 +
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
                                 broken_emp_t0 +
                                 j2has_dv_t0 +
                                 #rel_pov_t0 +
                                 health_t0 +
                                 srh_bin_t0 +
                                 ghq_case4_t0 +
                                 sf12mcs_dv_t0 +
                                 sf12pcs_dv_t0, 
                               data = not_pov_matchit_df,
                               method = "nearest", 
                               ratio=3,
                               distance = "glm",
                               estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

not_pov_matchit_df$weights_ps <- not_pov_matchit_mod$weights

summary(not_pov_matchit_df$weights_ps)

not_pov_matchit_df <- not_pov_matchit_df %>% filter(weights_ps!=0)

### create weighted data
not_pov_matchit_df_svy <- svydesign(ids = ~1,
                                    data = not_pov_matchit_df,
                                    weights = ~weights_ps)


### PS matched table one
not_pov_table_one_matchit <- svyCreateTableOne(vars = cov_vector3,
                                               strata = "exposure1",
                                               data = not_pov_matchit_df_svy,
                                               factorVars = c(catVars_vec),
                                               test = FALSE)

not_pov_table_one_matchit_sav <- print(not_pov_table_one_matchit, showAllLevels = TRUE, smd = TRUE,
                                       nonnormal = nonnorm_vec,
                                       factorVars = c(catVars_vec),
                                       formatOptions = list(big.mark = ","))


write.csv(not_pov_table_one_matchit_sav, "./output/cc/matched_descriptives/subgroups/rel_poverty/not_pov_table_one_matchit_sav.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(not_pov_table_one_matchit) > 0.1))

not_pov_table_one_matchit_smd <- data.frame(ExtractSmd(not_pov_table_one_matchit))
not_pov_table_one_matchit_smd <- not_pov_table_one_matchit_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(not_pov_table_one_matchit_smd, "./working_data/cc/subgroup/rel_poverty/not_pov_table_one_matchit_smd.csv")



################################################################################
#####                      IPTW using WeightIt package                     #####
################################################################################

#### rel_pov --------------------------------------------------------------------
### convert SF-12 outcomes to numeric to allow svyglm to work
rel_pov_weightit_df <- rel_pov_df  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed as in PS matching

### convert SF-12 outcomes to numeric to allow svyglm to work
rel_pov_weightit_df$sf12pcs_dv_t0 <- as.numeric(rel_pov_weightit_df$sf12pcs_dv_t0)
rel_pov_weightit_df$sf12mcs_dv_t0 <- as.numeric(rel_pov_weightit_df$sf12mcs_dv_t0)
rel_pov_weightit_df$sf12pcs_dv_t1 <- as.numeric(rel_pov_weightit_df$sf12pcs_dv_t1)
rel_pov_weightit_df$sf12mcs_dv_t1 <- as.numeric(rel_pov_weightit_df$sf12mcs_dv_t1)

rel_pov_weightit_df$srh_bin_t1 <- as.character(rel_pov_weightit_df$srh_bin_t1)
rel_pov_weightit_df$ghq_case4_t1 <- as.character(rel_pov_weightit_df$ghq_case4_t1)


start_time <- Sys.time()
rel_pov_weight_weightit <- weightit(exposure1 ~
                                      sex_dv_t0 +
                                      age_dv_t0 +
                                      non_white_t0 +
                                      #        marital_status_t0_married_civil_partnership +
                                      marital_status_t0_divorced_separated_widowed +
                                      marital_status_t0_single +
                                      dep_child_bin_t0 +
                                      degree_bin_t0 +
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
                                      broken_emp_t0 +
                                      j2has_dv_t0 +
                                      #rel_pov_t0 +
                                      health_t0 +
                                      srh_bin_t0 +
                                      ghq_case4_t0 +
                                      sf12mcs_dv_t0 +
                                      sf12pcs_dv_t0, 
                                    data = rel_pov_weightit_df, 
                                    stabilize = TRUE,
                                    estimand = "ATE",  
                                    method = "ps")  
rel_pov_weight_weightit
end_time <- Sys.time()
end_time - start_time

summary(rel_pov_weight_weightit)
summary(rel_pov_weight_weightit$weights)
density(rel_pov_weight_weightit$weights)

rel_pov_weightit_df$weightit_ipw <- rel_pov_weight_weightit$weights

rel_pov_weightit_df %>% group_by(exposure1) %>% summarise(n=sum(weightit_ipw))##

### created IPTW dataframe
rel_pov_weightit_df_svy <- svydesign(ids = ~1,
                                     data = rel_pov_weightit_df,
                                     weights = ~weightit_ipw)

### create IPTW table one
rel_pov_table_one_weightit <- svyCreateTableOne(vars = cov_vector3,
                                                strata = "exposure1",
                                                data = rel_pov_weightit_df_svy,
                                                factorVars = c(catVars_vec),
                                                test = FALSE)

rel_pov_table_one_weightit_sav <- print(rel_pov_table_one_weightit, showAllLevels = TRUE, smd = TRUE,
                                        nonnormal = nonnorm_vec,
                                        factorVars = c(catVars_vec),
                                        formatOptions = list(big.mark = ","))

write.csv(rel_pov_table_one_weightit_sav, "./output/cc/weighted_descriptives/subgroups/rel_poverty/rel_pov_table_one_weightit_sav.csv")

### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(rel_pov_table_one_weightit) > 0.1))
addmargins(table(ExtractSmd(rel_pov_table_one_weightit) > 0.2))

#### not_pov ----------------------------------------------------------------------
### convert SF-12 outcomes to numeric to allow svyglm to work
not_pov_weightit_df <- not_pov_df  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed as in PS matching

### convert SF-12 outcomes to numeric to allow svyglm to work
not_pov_weightit_df$sf12pcs_dv_t0 <- as.numeric(not_pov_weightit_df$sf12pcs_dv_t0)
not_pov_weightit_df$sf12mcs_dv_t0 <- as.numeric(not_pov_weightit_df$sf12mcs_dv_t0)
not_pov_weightit_df$sf12pcs_dv_t1 <- as.numeric(not_pov_weightit_df$sf12pcs_dv_t1)
not_pov_weightit_df$sf12mcs_dv_t1 <- as.numeric(not_pov_weightit_df$sf12mcs_dv_t1)

not_pov_weightit_df$srh_bin_t1 <- as.character(not_pov_weightit_df$srh_bin_t1)
not_pov_weightit_df$ghq_case4_t1 <- as.character(not_pov_weightit_df$ghq_case4_t1)


start_time <- Sys.time()
not_pov_weight_weightit <- weightit(exposure1 ~
                                      sex_dv_t0 +
                                      age_dv_t0 +
                                      non_white_t0 +
                                      #        marital_status_t0_married_civil_partnership +
                                      marital_status_t0_divorced_separated_widowed +
                                      marital_status_t0_single +
                                      dep_child_bin_t0 +
                                      degree_bin_t0 +
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
                                      broken_emp_t0 +
                                      j2has_dv_t0 +
                                      #rel_pov_t0 +
                                      health_t0 +
                                      srh_bin_t0 +
                                      ghq_case4_t0 +
                                      sf12mcs_dv_t0 +
                                      sf12pcs_dv_t0, 
                                    data = not_pov_weightit_df, 
                                    stabilize = TRUE,
                                    estimand = "ATE",  
                                    method = "ps")  
not_pov_weight_weightit
end_time <- Sys.time()
end_time - start_time

summary(not_pov_weight_weightit)
summary(not_pov_weight_weightit$weights)
density(not_pov_weight_weightit$weights)

not_pov_weightit_df$weightit_ipw <- not_pov_weight_weightit$weights

not_pov_weightit_df %>% group_by(exposure1) %>% summarise(n=sum(weightit_ipw))##

### created IPTW dataframe
not_pov_weightit_df_svy <- svydesign(ids = ~1,
                                     data = not_pov_weightit_df,
                                     weights = ~weightit_ipw)

### create IPTW table one
not_pov_table_one_weightit <- svyCreateTableOne(vars = cov_vector3,
                                                strata = "exposure1",
                                                data = not_pov_weightit_df_svy,
                                                factorVars = c(catVars_vec),
                                                test = FALSE)

not_pov_table_one_weightit_sav <- print(not_pov_table_one_weightit, showAllLevels = TRUE, smd = TRUE,
                                        nonnormal = nonnorm_vec,
                                        factorVars = c(catVars_vec),
                                        formatOptions = list(big.mark = ","))

write.csv(not_pov_table_one_weightit_sav, "./output/cc/weighted_descriptives/subgroups/rel_poverty/not_pov_table_one_weightit_sav.csv")

### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(not_pov_table_one_weightit) > 0.1))
addmargins(table(ExtractSmd(not_pov_table_one_weightit) > 0.2))


################################################################################
#####                               outcomes                              ######
################################################################################

# double robust MSM

#### rel_pov --------------------------------------------------------------------

### SF-12 PCS -----------------------
start_time <- Sys.time()
dr_iptw_pcs_rel_pov_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                    exposure1 +
                                    sex_dv_t0 +
                                    age_dv_t0 +
                                    age_dv_t1 +
                                      non_white_t0 +
                                      #        marital_status_t0_married_civil_partnership +
                                      marital_status_t0_divorced_separated_widowed +
                                      marital_status_t0_single +
                                      dep_child_bin_t0 +
                                      degree_bin_t0 +
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
                                      broken_emp_t0 +
                                      j2has_dv_t0 +
                                      #rel_pov_t0 +
                                      health_t0 +
                                    health_t1 +
                                    sf12pcs_dv_t0 +
                                    # interaction terms
                                    sex_dv_t0*age_dv_t0 +
                                    #sex_dv_t0*rel_pov_t0 +
                                    # age_dv_t0*rel_pov_t0 +
                                    (1|pidp),
                                  weights = weightit_ipw,
                                  data = rel_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

dr_iptw_pcs_rel_pov <- summary(dr_iptw_pcs_rel_pov_mod)

#fixef(dr_iptw_pcs_glmmTMB_mod)


dr_iptw_pcs_rel_pov_df <- tidy(dr_iptw_pcs_rel_pov_mod)

## confidence intervals
dr_iptw_pcs_rel_pov_df_ci <- data.frame(confint(dr_iptw_pcs_rel_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_pcs_rel_pov_df_ci$term <- rownames(dr_iptw_pcs_rel_pov_df_ci)

## join dfs together
dr_iptw_pcs_rel_pov_df <- dr_iptw_pcs_rel_pov_df %>% 
  left_join(dr_iptw_pcs_rel_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------
start_time <- Sys.time()
dr_iptw_mcs_rel_pov_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                    exposure1 +
                                    sex_dv_t0 +
                                    age_dv_t0 +
                                    age_dv_t1 +
                                      non_white_t0 +
                                      #        marital_status_t0_married_civil_partnership +
                                      marital_status_t0_divorced_separated_widowed +
                                      marital_status_t0_single +
                                      dep_child_bin_t0 +
                                      degree_bin_t0 +
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
                                      broken_emp_t0 +
                                      j2has_dv_t0 +
                                      #rel_pov_t0 +
                                      health_t0 +
                                    health_t1 +
                                    sf12mcs_dv_t0 +
                                    # interaction terms
                                    sex_dv_t0*age_dv_t0 +
                                    #sex_dv_t0*rel_pov_t0 +
                                    # age_dv_t0*rel_pov_t0 +
                                    (1|pidp),
                                    weights = weightit_ipw,
                                    data = rel_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_mcs_rel_pov_mod)

diagnose(dr_iptw_mcs_rel_pov_mod)

dr_iptw_mcs_rel_pov_df <- tidy(dr_iptw_mcs_rel_pov_mod)

## confidence intervals
dr_iptw_mcs_rel_pov_df_ci <- data.frame(confint(dr_iptw_mcs_rel_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_mcs_rel_pov_df_ci$term <- rownames(dr_iptw_mcs_rel_pov_df_ci)

## join dfs together
dr_iptw_mcs_rel_pov_df <- dr_iptw_mcs_rel_pov_df %>% 
  left_join(dr_iptw_mcs_rel_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

rel_pov_weightit_df$srh_bin_t1 <- factor(rel_pov_weightit_df$srh_bin_t1,
                                levels = c("excellent/very good/good", 
                                           "fair/poor"))


start_time <- Sys.time()
dr_iptw_srh_rel_pov_mod <- glmmTMB( srh_bin_t1 ~
                                    exposure1 +
                                    sex_dv_t0 +
                                    age_dv_t0 +
                                    age_dv_t1 +
                                      non_white_t0 +
                                      #        marital_status_t0_married_civil_partnership +
                                      marital_status_t0_divorced_separated_widowed +
                                      marital_status_t0_single +
                                      dep_child_bin_t0 +
                                      degree_bin_t0 +
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
                                      broken_emp_t0 +
                                      j2has_dv_t0 +
                                      #rel_pov_t0 +
                                      health_t0 +
                                    health_t1 +
                                    srh_bin_t0 +
                                    # interaction terms
                                    sex_dv_t0*age_dv_t0 +
                                    #sex_dv_t0*rel_pov_t0 +
                                    # age_dv_t0*rel_pov_t0 +
                                    (1|pidp),
                                  family = binomial(link="logit"),
                                  weights = weightit_ipw,
                                  data = rel_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_srh_rel_pov_mod)
diagnose(dr_iptw_srh_rel_pov_mod)

dr_iptw_srh_rel_pov_df <- tidy(dr_iptw_srh_rel_pov_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_srh_rel_pov_df_ci <- data.frame(confint(dr_iptw_srh_rel_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_srh_rel_pov_df_ci$term <- rownames(dr_iptw_srh_rel_pov_df_ci)

## join dfs together
dr_iptw_srh_rel_pov_df <- dr_iptw_srh_rel_pov_df %>% 
  left_join(dr_iptw_srh_rel_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------

rel_pov_weightit_df$ghq_case4_t1 <- factor(rel_pov_weightit_df$ghq_case4_t1,
                                  levels = c("0-3", "4 or more"))


start_time <- Sys.time()
dr_iptw_ghq_rel_pov_mod <- glmmTMB( ghq_case4_t1 ~
                                    exposure1 +
                                    sex_dv_t0 +
                                    age_dv_t0 +
                                    age_dv_t1 +
                                      non_white_t0 +
                                      #        marital_status_t0_married_civil_partnership +
                                      marital_status_t0_divorced_separated_widowed +
                                      marital_status_t0_single +
                                      dep_child_bin_t0 +
                                      degree_bin_t0 +
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
                                      broken_emp_t0 +
                                      j2has_dv_t0 +
                                      #rel_pov_t0 +
                                      health_t0 +
                                    health_t1 +
                                    ghq_case4_t0 +
                                    # interaction terms
                                    sex_dv_t0*age_dv_t0 +
                                    #sex_dv_t0*rel_pov_t0 +
                                    #age_dv_t0*rel_pov_t0 +
                                    (1|pidp),
                                  family = binomial(link="logit"),
                                  weights = weightit_ipw,
                                  data = rel_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_ghq_rel_pov_mod)
diagnose(dr_iptw_ghq_rel_pov_mod)

dr_iptw_ghq_rel_pov_df <- tidy(dr_iptw_ghq_rel_pov_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_ghq_rel_pov_df_ci <- data.frame(confint(dr_iptw_ghq_rel_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_ghq_rel_pov_df_ci$term <- rownames(dr_iptw_ghq_rel_pov_df_ci)

## join dfs together
dr_iptw_ghq_rel_pov_df <- dr_iptw_ghq_rel_pov_df %>% 
  left_join(dr_iptw_ghq_rel_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### create single summary df for iptw double-robust cc outcomes ------------

dr_iptw_rel_pov_df <- dr_iptw_pcs_rel_pov_df %>% 
  bind_rows(dr_iptw_mcs_rel_pov_df,
            dr_iptw_srh_rel_pov_df,
            dr_iptw_ghq_rel_pov_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, effect, est_type, estimate, std.error, p.value, lci, uci)

write.csv(dr_iptw_rel_pov_df, "./output/cc/weighted_outcomes/sub_groups/rel_poverty/dr_iptw_rel_pov_df.csv")


diagnose(dr_iptw_pcs_rel_pov_mod)
diagnose(dr_iptw_mcs_rel_pov_mod)
diagnose(dr_iptw_srh_rel_pov_mod)
diagnose(dr_iptw_ghq_rel_pov_mod)




#### not_pov ----------------------------------------------------------------------

### SF-12 PCS -----------------------
start_time <- Sys.time()
dr_iptw_pcs_not_pov_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                  exposure1 +
                                  sex_dv_t0 +
                                  age_dv_t0 +
                                  age_dv_t1 +
                                    non_white_t0 +
                                    #        marital_status_t0_married_civil_partnership +
                                    marital_status_t0_divorced_separated_widowed +
                                    marital_status_t0_single +
                                    dep_child_bin_t0 +
                                    degree_bin_t0 +
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
                                    broken_emp_t0 +
                                    j2has_dv_t0 +
                                    #rel_pov_t0 +
                                    health_t0 +
                                  health_t1 +
                                  sf12pcs_dv_t0 +
                                  # interaction terms
                                  sex_dv_t0*age_dv_t0 +
                                  #sex_dv_t0*rel_pov_t0 +
                                  #age_dv_t0*rel_pov_t0 +
                                  (1|pidp),
                                  weights = weightit_ipw,
                                  data = not_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

dr_iptw_pcs_not_pov <- summary(dr_iptw_pcs_not_pov_mod)

#fixef(dr_iptw_pcs_glmmTMB_mod)


dr_iptw_pcs_not_pov_df <- tidy(dr_iptw_pcs_not_pov_mod)

## confidence intervals
dr_iptw_pcs_not_pov_df_ci <- data.frame(confint(dr_iptw_pcs_not_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_pcs_not_pov_df_ci$term <- rownames(dr_iptw_pcs_not_pov_df_ci)

## join dfs together
dr_iptw_pcs_not_pov_df <- dr_iptw_pcs_not_pov_df %>% 
  left_join(dr_iptw_pcs_not_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------
start_time <- Sys.time()
dr_iptw_mcs_not_pov_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                  exposure1 +
                                  sex_dv_t0 +
                                  age_dv_t0 +
                                  age_dv_t1 +
                                    non_white_t0 +
                                    #        marital_status_t0_married_civil_partnership +
                                    marital_status_t0_divorced_separated_widowed +
                                    marital_status_t0_single +
                                    dep_child_bin_t0 +
                                    degree_bin_t0 +
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
                                    broken_emp_t0 +
                                    j2has_dv_t0 +
                                    #rel_pov_t0 +
                                    health_t0 +
                                  health_t1 +
                                  sf12mcs_dv_t0 +
                                  # interaction terms
                                  sex_dv_t0*age_dv_t0 +
                                  #sex_dv_t0*rel_pov_t0 +
                                  #age_dv_t0*rel_pov_t0 +
                                  (1|pidp),
                                  weights = weightit_ipw,
                                  data = not_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_mcs_not_pov_mod)

diagnose(dr_iptw_mcs_not_pov_mod)

dr_iptw_mcs_not_pov_df <- tidy(dr_iptw_mcs_not_pov_mod)

## confidence intervals
dr_iptw_mcs_not_pov_df_ci <- data.frame(confint(dr_iptw_mcs_not_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_mcs_not_pov_df_ci$term <- rownames(dr_iptw_mcs_not_pov_df_ci)

## join dfs together
dr_iptw_mcs_not_pov_df <- dr_iptw_mcs_not_pov_df %>% 
  left_join(dr_iptw_mcs_not_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

not_pov_weightit_df$srh_bin_t1 <- factor(not_pov_weightit_df$srh_bin_t1,
                                levels = c("excellent/very good/good", 
                                           "fair/poor"))


start_time <- Sys.time()
dr_iptw_srh_not_pov_mod <- glmmTMB( srh_bin_t1 ~
                                  exposure1 +
                                  sex_dv_t0 +
                                  age_dv_t0 +
                                  age_dv_t1 +
                                    non_white_t0 +
                                    #        marital_status_t0_married_civil_partnership +
                                    marital_status_t0_divorced_separated_widowed +
                                    marital_status_t0_single +
                                    dep_child_bin_t0 +
                                    degree_bin_t0 +
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
                                    broken_emp_t0 +
                                    j2has_dv_t0 +
                                    #rel_pov_t0 +
                                    health_t0 +
                                  health_t1 +
                                  srh_bin_t0 +
                                  # interaction terms
                                  sex_dv_t0*age_dv_t0 +
                                  #sex_dv_t0*rel_pov_t0 +
                                  #age_dv_t0*rel_pov_t0 +
                                  (1|pidp),
                                family = binomial(link="logit"),
                                weights = weightit_ipw,
                                data = not_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_srh_not_pov_mod)
diagnose(dr_iptw_srh_not_pov_mod)

dr_iptw_srh_not_pov_df <- tidy(dr_iptw_srh_not_pov_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_srh_not_pov_df_ci <- data.frame(confint(dr_iptw_srh_not_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_srh_not_pov_df_ci$term <- rownames(dr_iptw_srh_not_pov_df_ci)

## join dfs together
dr_iptw_srh_not_pov_df <- dr_iptw_srh_not_pov_df %>% 
  left_join(dr_iptw_srh_not_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------

not_pov_weightit_df$ghq_case4_t1 <- factor(not_pov_weightit_df$ghq_case4_t1,
                                  levels = c("0-3", "4 or more"))


start_time <- Sys.time()
dr_iptw_ghq_not_pov_mod <- glmmTMB( ghq_case4_t1 ~
                                  exposure1 +
                                  sex_dv_t0 +
                                  age_dv_t0 +
                                  age_dv_t1 +
                                    non_white_t0 +
                                    #        marital_status_t0_married_civil_partnership +
                                    marital_status_t0_divorced_separated_widowed +
                                    marital_status_t0_single +
                                    dep_child_bin_t0 +
                                    degree_bin_t0 +
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
                                    broken_emp_t0 +
                                    j2has_dv_t0 +
                                    #rel_pov_t0 +
                                    health_t0 +
                                  health_t1 +
                                  ghq_case4_t0 +
                                  # interaction terms
                                  sex_dv_t0*age_dv_t0 +
                                  #sex_dv_t0*rel_pov_t0 +
                                  #age_dv_t0*rel_pov_t0 +
                                  (1|pidp),
                                family = binomial(link="logit"),
                                weights = weightit_ipw,
                                data = not_pov_weightit_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_ghq_not_pov_mod)
diagnose(dr_iptw_ghq_not_pov_mod)

dr_iptw_ghq_not_pov_df <- tidy(dr_iptw_ghq_not_pov_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_ghq_not_pov_df_ci <- data.frame(confint(dr_iptw_ghq_not_pov_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_ghq_not_pov_df_ci$term <- rownames(dr_iptw_ghq_not_pov_df_ci)

## join dfs together
dr_iptw_ghq_not_pov_df <- dr_iptw_ghq_not_pov_df %>% 
  left_join(dr_iptw_ghq_not_pov_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### create single summary df for iptw double-robust cc outcomes ------------

dr_iptw_not_pov_df <- dr_iptw_pcs_not_pov_df %>% 
  bind_rows(dr_iptw_mcs_not_pov_df,
            dr_iptw_srh_not_pov_df,
            dr_iptw_ghq_not_pov_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, effect, est_type, estimate, std.error, p.value, lci, uci)

write.csv(dr_iptw_not_pov_df, "./output/cc/weighted_outcomes/sub_groups/rel_poverty/dr_iptw_not_pov_df.csv")


diagnose(dr_iptw_pcs_not_pov_mod)
diagnose(dr_iptw_mcs_not_pov_mod)
diagnose(dr_iptw_srh_not_pov_mod)
diagnose(dr_iptw_ghq_not_pov_mod)
