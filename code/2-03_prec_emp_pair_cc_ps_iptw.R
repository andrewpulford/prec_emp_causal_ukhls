################################################################################

# Precarious employment and health - Understanding Society
# 3-04 - Paired unweighted complete case analysis propensity scores for risk of 
# job loss 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) produces unweighted risk of job loss propensity scores for complete case 
#     paired data


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
library(WeightIt) # for inverse probability weighting
library(tableone) # for creating table one
library(survey) # for PS weighting
library(reshape2) # Reorganizing data
library(broom) # for tidying regression outputs into df format
library(lme4) # for multi-level modelling
library(glmmTMB) # for multi-level modelling (faster than lme4)
library(numDeriv) # for gradient checks

################################################################################
#####                         load and prepare data                        #####
################################################################################

####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")

## numeric non-normal vars
nonnorm_vec <- (c("age_dv_t0", "sf12mcs_dv_t0", "sf12pcs_dv_t0"))

####load eligible cases --------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/cc/pair_cc_analytic.rds")

### convert binary outcome vars to factors to allow svyglm to work
pair_cc_analytic$srh_bin_t1 <- factor(pair_cc_analytic$srh_bin_t1)
pair_cc_analytic$ghq_case4_t1 <- factor(pair_cc_analytic$ghq_case4_t1)

### convert SF-12 outcomes to numeric to allow svyglm to work
pair_cc_analytic$sf12pcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t1)
pair_cc_analytic$sf12mcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t1)



#### keep only variables required for propensity score (and other analysis) ----
pair_cc_ps <- pair_cc_analytic %>% 
  dplyr::select(pidp, all_of(c(cov_vector, cov_vector2, outcome_vector))) #%>% 
#  dplyr::select(-c(psu, strata, wt_name, wt_value))

#### convert pidp to factor 
pair_cc_ps$pidp <- factor(pair_cc_ps$pidp)

#### convert scale vars to numeric ---------------------------------------------
pair_cc_ps <- pair_cc_ps %>% 
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
pair_cc_ps[sapply(pair_cc_ps, is.character)] <- lapply(pair_cc_ps[sapply(pair_cc_ps, is.character)], as.factor)

#### revel exposure vars so that unexposed is lower ref category ---------------
pair_cc_ps$exposure1 <- factor(pair_cc_ps$exposure1,
                                     levels = rev(levels(pair_cc_ps$exposure1)))

pair_cc_ps$exposure2 <- factor(pair_cc_ps$exposure2,
                                     levels = rev(levels(pair_cc_ps$exposure2)))

#### create dummy variables for categorical vars with>2 cats -------------------

pair_cc_ps <- pair_cc_ps %>% 
  dummy_cols(select_columns = c("sex_dv_t0",
                                "non_white_t0",
                                "marital_status_t0",
                                "degree_bin_t0",
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
                                "ghq_case4_t0"
                                ))

pair_cc_ps <- pair_cc_ps %>% janitor::clean_names()

### use this if full df is too big
#test_df <- pair_cc_ps %>% 
#  slice_sample(prop = 1)

#### load unmatched SMD df
table_one_unmatched_smd <- read.csv("./working_data/cc/table_one_unmatched_smd.csv") %>% 
  dplyr::select(c(var, smd, imbalance_flag, matched))


################################################################################
#####                               functions                              #####
################################################################################

#### propensity score model - basic --------------------------------------------
ps_model <- function(data = pair_cc_ps, outcome){
  glm(outcome ~
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
          rel_pov_t0 +
          health_t0 +
          srh_bin_t0 +
          ghq_case4_t0 +
          sf12mcs_dv_t0 +
          sf12pcs_dv_t0,
        family = binomial(link="logit"),
        data = data)
  
  
}

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
#####           Propensity score matching using MatchIt package            #####
################################################################################

matchit_df <- pair_cc_ps  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed to allow 3:1 ratio matching

### convert SF-12 outcomes to numeric to allow svyglm to work
matchit_df$sf12pcs_dv_t0 <- as.numeric(matchit_df$sf12pcs_dv_t0)
matchit_df$sf12mcs_dv_t0 <- as.numeric(matchit_df$sf12mcs_dv_t0)
matchit_df$sf12pcs_dv_t1 <- as.numeric(matchit_df$sf12pcs_dv_t1)
matchit_df$sf12mcs_dv_t1 <- as.numeric(matchit_df$sf12mcs_dv_t1)

matchit_df$srh_bin_t1 <- as.character(matchit_df$srh_bin_t1)
matchit_df$ghq_case4_t1 <- as.character(matchit_df$ghq_case4_t1)


start_time <- Sys.time()
matchit_mod <- matchit(exp1_bin ~
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
                  rel_pov_t0 +
                  health_t0 +
                  srh_bin_t0 +
                  ghq_case4_t0 +
                  sf12mcs_dv_t0 +
                  sf12pcs_dv_t0, 
                data = matchit_df,
                method = "nearest", 
                ratio=3,
                distance = "glm",
                estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

matchit_df$weights_ps <- matchit_mod$weights

### create weighted data
matchit_df_svy <- svydesign(ids = ~1,
                            data = matchit_df,
                            weights = ~weights_ps)


### PS matched table one
table_one_matchit <- svyCreateTableOne(vars = cov_vector3,
                                        strata = "exposure1",
                                        data = matchit_df_svy,
                                       factorVars = c(catVars_vec),
                                       test = FALSE)

table_one_matchit_sav <- print(table_one_matchit, showAllLevels = TRUE, smd = TRUE,
                                nonnormal = nonnorm_vec,
                               factorVars = c(catVars_vec),
                               formatOptions = list(big.mark = ","))

write.csv(table_one_matchit_sav, "./output/temp_output/table_one_MatchIt_sav.csv")

### shorter version of table one for paper ------------
### matched table one
table_one_matchit2 <- svyCreateTableOne(vars = cov_vector,
                                       strata = "exposure1",
                                       data = matchit_df_svy,
                                       factorVars = c(catVars_vec),
                                       test = FALSE)

table_one_matchit2_sav <- print(table_one_matchit2, showAllLevels = TRUE, smd = FALSE,
                               nonnormal = nonnorm_vec,
                               factorVars = c(catVars_vec),
                               formatOptions = list(big.mark = ","))

write.csv(table_one_matchit2_sav, "./output/cc/matched_descriptives/table_one_IPTW_paper.csv")

### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(table_one_matchit) > 0.1))
addmargins(table(ExtractSmd(table_one_matchit) > 0.2))












################################################################################
#####                     IPTW using WeightIt package                      #####
################################################################################

weightit_df <- pair_cc_ps  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed as in PS matching

### convert SF-12 outcomes to numeric to allow svyglm to work
weightit_df$sf12pcs_dv_t0 <- as.numeric(weightit_df$sf12pcs_dv_t0)
weightit_df$sf12mcs_dv_t0 <- as.numeric(weightit_df$sf12mcs_dv_t0)
weightit_df$sf12pcs_dv_t1 <- as.numeric(weightit_df$sf12pcs_dv_t1)
weightit_df$sf12mcs_dv_t1 <- as.numeric(weightit_df$sf12mcs_dv_t1)

weightit_df$srh_bin_t1 <- as.character(weightit_df$srh_bin_t1)
weightit_df$ghq_case4_t1 <- as.character(weightit_df$ghq_case4_t1)


start_time <- Sys.time()
weight_weightit <- weightit(exposure1 ~
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
                              rel_pov_t0 +
                              health_t0 +
                              srh_bin_t0 +
                              ghq_case4_t0 +
                              sf12mcs_dv_t0 +
                              sf12pcs_dv_t0 ,
                                     data = weightit_df, 
                                     stabilize = TRUE,
                                     estimand = "ATE",  
                                     method = "ps")  
weight_weightit
end_time <- Sys.time()
end_time - start_time

summary(weight_weightit)
summary(weight_weightit$weights)
density(weight_weightit$weights)

weightit_df$weightit_ipw <- weight_weightit$weights

weightit_df %>% group_by(exposure1) %>% summarise(n=sum(weightit_ipw))##

### created IPTW dataframe
weightit_df_svy <- svydesign(ids = ~1,
                             data = weightit_df,
                             weights = ~weightit_ipw)

### create IPTW table one
table_one_weightit <- svyCreateTableOne(vars = cov_vector3,
                                   strata = "exposure1",
                                   data = weightit_df_svy,
                                   factorVars = c(catVars_vec),
                                   test = FALSE)

table_one_weightit_sav <- print(table_one_weightit, showAllLevels = TRUE, smd = TRUE,
                           nonnormal = nonnorm_vec,
                           factorVars = c(catVars_vec),
                           formatOptions = list(big.mark = ","))

write.csv(table_one_weightit_sav, "./output/cc/weighted_descriptives/table_one_weightit_sav.csv")

### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(table_one_weightit) > 0.1))
addmargins(table(ExtractSmd(table_one_weightit) > 0.2))


###### write working files -------------------------------------------------------
write_rds(matchit_df, "working_data/cc/matchit_df.rds") # propesntiy score matched analytic df (MatchIt)
write_rds(weightit_df, "working_data/cc/weightit_df.rds") # IPTW analytic df (WeightIt)


################################################################################
#####                  assess balance pre/post-matching                   ######
################################################################################

#### create df for matched smd values ------------------------------------------
table_one_matchit_smd <- data.frame(ExtractSmd(table_one_matchit))
table_one_matchit_smd <- table_one_matchit_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "Propensity score matched")


#### create df for IPW weighted smd values -----------------------------------------
### WeightIt --------
table_one_WeightIt_smd <- data.frame(ExtractSmd(table_one_weightit))
table_one_WeightIt_smd <- table_one_WeightIt_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "Inverse probability of treatment weighted")

#### bind to unmatched smd df --------------------------------------------------
assess_matching_balance <- table_one_unmatched_smd %>% 
  bind_rows(table_one_matchit_smd, 
            table_one_WeightIt_smd) %>% 
  mutate(smd_flag = ifelse(smd>0.2,1,0)) %>% 
  rename("method" = "matched")

### plot
balance_plot <- assess_matching_balance %>% 
  ggplot(aes(x=smd, y=var, col=method, shape=method)) +
  geom_point() +
#  geom_jitter(height = 0.01, width = 0.01) +
  geom_vline(xintercept = 0.1, linetype = "dotted") +
  geom_vline(xintercept = 0.2, linetype = "dashed") +
#  geom_text(aes(label=smd_flag), size = 3, position = "dodge") +
  theme_bw() +
  scale_color_manual(values = c("blue","red", "green"))

## id vars over thresholds
assess_matching_balance %>% filter(smd_flag==1) # SMD>0.2 
assess_matching_balance %>% filter(imbalance_flag=="SMD>0.1" & smd_flag==0) # 0.1>SMD<0.2 

tiff("./output/cc/balance_plot.tiff", width = 960)
balance_plot
dev.off()

#### number of individuals -----------------------------------------------------
#iptw_df %>% dplyr::select(pidp) %>% unique() %>% nrow()

