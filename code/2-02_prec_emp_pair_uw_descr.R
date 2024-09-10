################################################################################

# Precarious employment and health - Understanding Society
# 2-02 - Paired unweighted descriptives 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) produces unweighted descriptive analysis for complete case paired data


################################################################################

## remove any existing objects from global environment
rm(list=ls()) 

## turn off scientific notation
options(scipen = 999)

################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
library(Hmisc) # histogram plotting across cols
library(tableone) # for creating table one
library(fastDummies)
library(ggsankeyfier)
################################################################################
#####                         load and prepare data                        #####
################################################################################

####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")

####load analytic df -----------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/cc/pair_cc_analytic.rds")

# check no NAs
sapply(pair_cc_analytic, function(x) sum(is.na(x)))

#pair_cc_analytic$fimnnet_dv_t0 <- as.numeric(as.character(pair_cc_analytic$fimnnet_dv_t0))


#pair_cc_analytic$srh_dv_t0 <- factor(pair_cc_analytic$srh_dv_t0, 
#                                     levels = c("excellent", "very good", "good",
#                                                "fair", "poor"))

pair_cc_analytic$sf12pcs_dv_t0 <- as.character(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t0)

pair_cc_analytic$sf12mcs_dv_t0 <- as.character(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t0)

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
                                "ghq_case4_t0"
  ))

pair_cc_analytic <- pair_cc_analytic %>% janitor::clean_names()


#####----------------------------------------------------------------------#####
#####                          Save sample chars data                      #####
#####----------------------------------------------------------------------#####

## as dataframe
#write_rds(sample_chars, "./working_data/sample_chars.rds")
#
#
#
##### job loss between t0-t1
#
#emp_status_t1 <- pair_cc_analytic %>%
#  group_by(exposure1)  %>%  
#  summarise(n=n()) %>%  
#  mutate(est = n/sum(n)*100) %>% 
#  mutate(var="Employment status t1") %>% 
#  rename("measure"= "exposure1") %>% 
#  dplyr::select(var, measure, n, est)
#
## if unemployed at t0 or if unemp spell between t0-1
#lost_job <- pair_cc_analytic %>% 
#  group_by(exposure2) %>% 
#  summarise(n=n()) %>%  
#  mutate(est = n/sum(n)*100) %>% 
#  mutate(var="Job loss between t0 and t1") %>% 
#  rename("measure"= "exposure2") %>% 
#  dplyr::select(var, measure, n, est)
#
#
################################################################################
#####        Table 1: participant characteristics by exposure group        #####
################################################################################

pair_cc_analytic <- droplevels(pair_cc_analytic)

#pair_cc_analytic <- pair_cc_analytic %>% 
#  mutate(across(.cols = everything(), 
#                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 

pair_cc_analytic <- pair_cc_analytic %>% 
  dplyr::select(all_of(cov_vector3), exposure1, exposure2)

## histogram for each numeric variable
# temp df with only numeric vars
temp <- pair_cc_analytic %>% dplyr::select(c(age_dv_t0,
                                             sf12pcs_dv_t0,
                                             sf12mcs_dv_t0))

# histogram
#hist.data.frame(temp)

# normality test for age
temp2 <- sample (temp$age_dv_t0, size=5000)
shapiro.test(temp2)

# winsorized version of income to deal with extreme values
#hist(DescTools::Winsorize(pair_cc_analytic$fimnnet_dv_t0))

## vector for numeric vars
nonnorm_vec <- colnames(temp)


## unemployed at T1
table_one <- tableone::CreateTableOne(vars = cov_vector3, strata = "exposure1",
                            data = pair_cc_analytic,
                            factorVars = c(catVars_vec),
                            test = FALSE)

table_one_sav <- print(table_one, showAllLevels = TRUE, smd = TRUE,
                       nonnormal = nonnorm_vec,
                       factorVars = c(catVars_vec),
                       formatOptions = list(big.mark = ","))

# Count covariates with important imbalance
table_one_smd <- data.frame(ExtractSmd(table_one))
table_one_smd <- table_one_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_smd, "./working_data/cc/table_one_unmatched_smd.csv")

addmargins(table(ExtractSmd(table_one) > 0.1))

# plot unmatched smd's
table_one_smd %>% 
  ggplot(aes(x=smd, y=var, col=imbalance_flag)) +
  geom_point() +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("blue","red"))

## job loss between t0 and t1
table_one_alt <- CreateTableOne(vars = cov_vector3, strata = "exposure2", 
                                data = pair_cc_analytic,
                                factorVars = c(catVars_vec))

# Count covariates with important imbalance
table_one_alt_smd <- data.frame(ExtractSmd(table_one_alt))
table_one_alt_smd <- table_one_alt_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_alt_smd, "./output/cc/unmatched_descriptives/table_one_alt_unmatched_smd.csv")

# plot unmatched smd's
table_one_alt_smd %>% 
  ggplot(aes(x=smd, y=var, col=imbalance_flag)) +
  geom_point() +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("blue","red"))

addmargins(table(ExtractSmd(table_one_alt) > 0.1))

table_one_alt_sav <- print(table_one_alt, showAllLevels = TRUE, smd = TRUE,
                           nonnormal = nonnorm_vec,
                           factorVars = c(catVars_vec),
                           formatOptions = list(big.mark = ","))

### save tables

write.csv(table_one_sav, "./output/cc/unmatched_descriptives/table_one_unmatched.csv")
write.csv(table_one_alt_sav, "./output/cc/unmatched_descriptives/table_one_alt_unmatched.csv")

################################################################################
#####                       unweighted outcome models                      #####
################################################################################

library(glmmTMB) # for multi-level modelling (faster than lme4)

#### load IPTW df --------------------------------------------------------------
# use without weights
iptw_df <- readRDS("working_data/cc/weightit_df.rds") # analytic df with IPTW from WeightIt package

iptw_df <- iptw_df %>% 
  dummy_cols(select_columns = c("marital_status_t1",
                                "health_t1"
  ))

iptw_df <- iptw_df %>% janitor::clean_names()


### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
iptw_df$srh_bin_t0 <- factor(iptw_df$srh_bin_t0,
                             levels = c("excellent/very good", 
                                        "good/fair/poor"))
iptw_df$srh_bin_t1 <- factor(iptw_df$srh_bin_t1,
                             levels = c("excellent/very good", 
                                        "good/fair/poor"))


iptw_df$ghq_case4_t0 <- factor(iptw_df$ghq_case4_t0,
                               levels = c("0-3", "4 or more"))
iptw_df$ghq_case4_t1 <- factor(iptw_df$ghq_case4_t1,
                               levels = c("0-3", "4 or more"))

iptw_df <- iptw_df %>% 
  mutate(ghq4_outcome_bin = ifelse(ghq_case4_t1=="0-3",0,
                                   ifelse(ghq_case4_t1=="4 or more",1,
                                          "CHECK")))

iptw_df$exposure1 <- factor(iptw_df$exposure1,
                            levels = c("unexposed",
                                       "exposed (employed at t1)"))
iptw_df$exposure2 <- factor(iptw_df$exposure2,
                            levels = c("unexposed",
                                       "exposed (no job loss between t0 and t1"))
# error with exposure 2 to correct at some point

### convert SF-12 outcomes to numeric to allow svyglm to work
iptw_df$sf12pcs_dv_t0 <- as.numeric(iptw_df$sf12pcs_dv_t0)
iptw_df$sf12mcs_dv_t0 <- as.numeric(iptw_df$sf12mcs_dv_t0)
iptw_df$sf12pcs_dv_t1 <- as.numeric(iptw_df$sf12pcs_dv_t1)
iptw_df$sf12mcs_dv_t1 <- as.numeric(iptw_df$sf12mcs_dv_t1)

### ghq-12 change var for checking findings
iptw_df <- iptw_df %>% 
  mutate(ghq_change = paste0(ghq_case4_t0," to ", ghq_case4_t1))

table(iptw_df$ghq_change)

#### GHQ-12 --------------------------------------------------------------------
start_time <- Sys.time()
dr_ghq_glmmTMB_mod <- glmmTMB( ghq_case4_t1 ~
                                         exposure1 +
                                         ghq_case4_t0 +
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
                                         health_t1 +
                                         # interaction terms
                                         sex_dv_t0*age_dv_t0 +
                                         sex_dv_t0*rel_pov_t0 +
                                         age_dv_t0*rel_pov_t0 +
                                         (1|pidp),
                                       family = binomial(link="logit"),
                                       data = iptw_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_ghq_glmmTMB_mod)


################################################################################
##### scrapbook #####
################################################################################

##### sankeys GHQ-12 -----------------------------------------------------------

#### exposed group
temp1 <- pair_cc_analytic %>% 
  group_by(as.numeric(scghq2_dv_t0),as.numeric(scghq2_dv_t1),exposure1) %>% 
  summarise(n=n()) %>% 
  rename("ghq2_score_t0" = `as.numeric(scghq2_dv_t0)`,
         "ghq2_score_t1" = `as.numeric(scghq2_dv_t1)`) %>% 
  ungroup() %>% 
  group_by(exposure1) %>% 
  mutate(d=sum(n),
         p_t0 = n/d*100,
         code = paste0(ghq2_score_t0,"[",format(round(p_t0, 2), nsmall = 2),"]",ghq2_score_t1,"a"))


temp2 <- temp1 %>% filter(exposure1=="exposed (employed at t1)")

write.csv(temp2, "./output/scrapbook/ghq_sankey_treated.csv")

temp3 <- temp1 %>% filter(exposure1=="unexposed")

write.csv(temp3, "./output/scrapbook/ghq_sankey_control.csv")

###### basic regression model fo GHQ-12 ----------------------------------------

#### GHQ-12 --------------------------------------------------------------------
start_time <- Sys.time()
dr_ghq_glmmTMB_mod <- glmmTMB( ghq_case4_t1 ~
                                 exposure1,# +
#                                 ghq_case4_t0 +
#                                 sex_dv_t0 +
#                                 age_dv_t0 +
#                                 age_dv_t1 +
#                                 non_white_t0 +
#                                 #        marital_status_t0_married_civil_partnership +
#                                 marital_status_t0_divorced_separated_widowed +
#                                 marital_status_t0_single +
#                                 #        marital_status_t1_married_civil_partnership +
#                                 marital_status_t1_divorced_separated_widowed +
#                                 marital_status_t1_single +
#                                 dep_child_bin_t0 +
#                                 degree_bin_t0 +
#                                 #  gor_dv_t0_east_midlands +
#                                 gor_dv_t0_east_of_england +
#                                 gor_dv_t0_london +
#                                 gor_dv_t0_north_east +
#                                 gor_dv_t0_north_west +
#                                 gor_dv_t0_northern_ireland +
#                                 gor_dv_t0_scotland +
#                                 gor_dv_t0_south_east +
#                                 gor_dv_t0_south_west +
#                                 gor_dv_t0_wales +
#                                 gor_dv_t0_west_midlands +
#                                 gor_dv_t0_yorkshire_and_the_humber +
#                                 #  sic2007_section_lab_t0_accommodation_and_food_service_activities +
#                                 sic2007_section_lab_t0_administrative_and_support_service_activities +
#                                 sic2007_section_lab_t0_construction +
#                                 sic2007_section_lab_t0_education +
#                                 sic2007_section_lab_t0_human_health_and_social_work_activities +
#                                 sic2007_section_lab_t0_manufacturing +
#                                 sic2007_section_lab_t0_other_industry +
#                                 sic2007_section_lab_t0_professional_scientific_and_technical_activities +
#                                 sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security +
#                                 sic2007_section_lab_t0_transportation_and_storage +
#                                 sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles +
#                                 #  soc2000_major_group_title_t0_administrative_and_secretarial_occupations +
#                                 soc2000_major_group_title_t0_associate_professional_and_technical_occupations +
#                                 soc2000_major_group_title_t0_elementary_occupations +
#                                 soc2000_major_group_title_t0_managers_and_senior_officials +
#                                 soc2000_major_group_title_t0_personal_service_occupations +
#                                 soc2000_major_group_title_t0_process_plant_and_machine_operatives +
#                                 soc2000_major_group_title_t0_sales_and_customer_service_occupations +
#                                 soc2000_major_group_title_t0_science_and_technology_professionals +
#                                 soc2000_major_group_title_t0_skilled_trades_occupations +
#                                 jbft_dv_t0 +
#                                 small_firm_t0 +
#                                 emp_contract_t0 +
#                                 broken_emp_t0 +
#                                 j2has_dv_t0 +
#                                 rel_pov_t0 +
#                                 health_t0 +
#                                 health_t1 +
#                                 # interaction terms
#                                 sex_dv_t0*age_dv_t0 +
#                                 sex_dv_t0*rel_pov_t0 +
#                                 age_dv_t0*rel_pov_t0 +
#                                 (1|pidp),
                               family = binomial(link="logit"),
                               data = iptw_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_ghq_glmmTMB_mod)


###### linear prob model

lm(ghq_case4_t1 ~ exposure1,
  data=iptw_df)
