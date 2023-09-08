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
library(Matching) # PS matching
library(tableone) # for creating table one
library(survey) # for PS weighting
library(reshape2) # Reorganizing data
library(broom) # for tidying regression outputs into df format
library(lme4) # for multi-level modelling
library(glmmTMB) # for multi-level modelling (faster than lme4)
library(numDeriv) # for gradient checks
library(ipw) # for inverse probability weighting

################################################################################
#####                         load and prepare data                        #####
################################################################################

####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")

## numeric non-normal vars
nonnorm_vec <- (c("age_dv_t0", "sf12mcs_dv_t0", "sf12pcs_dv_t0"))

####load eligible cases --------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")

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

pair_cc_ps <- pair_cc_ps %>% janitor::clean_names()

#### remove any NAs not coded as a missing category ----------------------------
##  check NAs
#sapply(pair_cc_ps, function(x) sum(is.na(x)))
## remove
#pair_cc_ps <- na.omit(pair_cc_ps)
## check again
#sapply(pair_cc_ps, function(x) sum(is.na(x)))


### use this if full df is too big
#test_df <- pair_cc_ps %>% 
#  slice_sample(prop = 1)

#### load unmatched SMD df
table_one_unmatched_smd <- read.csv("./working_data/table_one_unmatched_smd.csv") %>% 
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
#####                      IPTW using MatchIt package                      #####
################################################################################

matchit_df <- pair_cc_ps  %>%  
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           1,0))

### convert SF-12 outcomes to numeric to allow svyglm to work
pair_cc_analytic$sf12pcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t1)
pair_cc_analytic$sf12mcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t1)

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
                data = matchit_df,
                method = "quick", # Generalized Full Matching
                distance = "glm",
                estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

#test_summary <- summary(test)
#write.csv(test_summary, "./output/temp_output/test_sumary.csv")

matchit_df$weights_ps <- matchit_mod$weights
#matchit_df$weights <- unname(matchit_df$weights)

### create weighted data
matchit_df_svy <- svydesign(ids = ~1,
                            data = matchit_df,
                            weights = ~weights_ps)


### weighted table one
table_one_matchit <- svyCreateTableOne(vars = cov_vector3,
                                        strata = "exposure1",
                                        data = matchit_df_svy,
                                        test = FALSE)

table_one_matchit_sav <- print(table_one_matchit, showAllLevels = TRUE, smd = TRUE,
                                nonnormal = nonnorm_vec,
                                formatOptions = list(big.mark = ","))

write.csv(table_one_matchit_sav, "./output/temp_output/table_one_MatchItQuick_sav.csv")

### shorter version of table one for paper ------------
### weighted table one
table_one_matchit2 <- svyCreateTableOne(vars = cov_vector,
                                       strata = "exposure1",
                                       data = matchit_df_svy,
                                       test = FALSE)

table_one_matchit2_sav <- print(table_one_matchit2, showAllLevels = TRUE, smd = FALSE,
                               nonnormal = nonnorm_vec,
                               formatOptions = list(big.mark = ","))

write.csv(table_one_matchit2_sav, "./output/weighted_descriptives/table_one_IPTW_paper.csv")

### count covariates with an important imbalance (>0.1 or >0.2)
addmargins(table(ExtractSmd(table_one_matchit) > 0.1))
addmargins(table(ExtractSmd(table_one_matchit) > 0.2))


################################################################################
#####                     Unemployment at t1 PS model                      #####
################################################################################

### call the function (and benchmark time - takes a while to run for MLM 3 hours)
## not MLM
start_time <- Sys.time()
ps_mod_exp1_noMLM <- ps_model(data = pair_cc_ps, outcome = pair_cc_ps$exposure1)
end_time <- Sys.time()
end_time - start_time

summary(ps_mod_exp1_noMLM)

## MLM
start_time <- Sys.time()
ps_mod_exp1 <- ps_model_mlm(data = pair_cc_ps, outcome = pair_cc_ps$exposure1)
end_time <- Sys.time()
end_time - start_time

### summary of model
summary(ps_mod_exp1)

### checks
#
## rescale vars
#numcols <- grep("^c\\.",names(pair_cc_ps))
#pair_cc_ps_scaled <- pair_cc_ps
#pair_cc_ps_scaled[,numcols] <- scale(pair_cc_ps_scaled[,numcols])
#ps_mod_exp1 <- update(ps_mod_exp1,data=pair_cc_ps_scaled)#
#
## Check singularity
#tt <- getME(ps_mod_exp1,"theta")
#ll <- getME(ps_mod_exp1,"lower")
#min(tt[ll==0])
#
## Double-checking gradient calculations
#derivs1 <- ps_mod_exp1@optinfo$derivs
#check_grad1 <- with(derivs1,solve(Hessian,gradient))
#max(abs(check_grad1))
#
# try with numDeriv
#dd <- update(ps_mod_exp1,devFunOnly=TRUE)
#pars <- unlist(getME(ps_mod_exp1,c("theta","fixef")))
#grad2 <- grad(dd,pars)
#hess2 <- hessian(dd,pars)
#check_grad2 <- solve(hess2,grad2)
#max(pmin(abs(check_grad2),abs(grad2)))

## try update with max iterations
#ss <- getME(ps_mod_exp1,c("theta","fixef"))
#start_time <- Sys.time()
#ps_mod_exp1_update <- update(ps_mod_exp1,start=ss,
#                             control=glmerControl(optimizer="bobyqa",
#                                                  optCtrl=list(maxfun=2e5)))
#end_time <- Sys.time()
#end_time - start_time



## check for 1 or 0 predicted probabilities
pair_cc_ps$y_pred <- predict(ps_mod_exp1, pair_cc_ps, type="response")

summary(pair_cc_ps$y_pred)

### predicted probability of being assigned to exposed group
pair_cc_ps$ps_exp1 <- predict(ps_mod_exp1, type = "response")
summary(pair_cc_ps$ps_exp1)

### predicted probability of being assigned to unexposed group
pair_cc_ps$ps_noexp1 <- 1-pair_cc_ps$ps_exp1
summary(pair_cc_ps$ps_noexp1)

### predicted probability of being assigned to actual exposure status
pair_cc_ps$ps_assign <- NA
pair_cc_ps$ps_assign[pair_cc_ps$exposure1=="exposed (employed at t1)"] <- pair_cc_ps$ps_exp1[pair_cc_ps$exposure1=="exposed (employed at t1)"]
pair_cc_ps$ps_assign[pair_cc_ps$exposure1=="unexposed"] <- pair_cc_ps$ps_noexp1[pair_cc_ps$exposure1=="unexposed"]

### smaller of ps_exp1 and ps_noexp1 for matchnig weight
pair_cc_ps$ps_min <- pmin(pair_cc_ps$ps_exp1, pair_cc_ps$ps_noexp1)
summary(pair_cc_ps$ps_min)

################################################################################
#####                       propensity score matching                      #####
################################################################################

#### match propensity scores using Matching package ----------------------------
list_match <- Match(Tr = (pair_cc_ps$exposure1=="exposed (employed at t1)"),
# logit of PS/1-PS
X = log(pair_cc_ps$ps_exp1/pair_cc_ps$ps_noexp1),
## 1:1 matching ratio
M = 1,
## caliper = 0.2 * SD(logit(PS))
caliper  = 0.2,
replace  = FALSE,
ties     = TRUE,
version  = "fast")

#### extract matched data ------------------------------------------------------
df_matched <- pair_cc_ps[unlist(list_match[c("index.treated","index.control")]), ]

#### matched Table One with SMD ------------------------------------------------
table_one_matched <- CreateTableOne(vars = cov_vector3,
                                    strata = "exposure1",
                                    data = df_matched,
                                    test = FALSE)

table_one_matched_sav <- print(table_one_matched,  smd = TRUE)

write.csv(table_one_matched_sav, "./output/matched_descriptives/table_one_matched.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(table_one_matched) > 0.1))


################################################################################
#####                   propensity score matching weight                   #####
################################################################################

### create matching weight
pair_cc_ps$mw <- pair_cc_ps$ps_min/pair_cc_ps$ps_assign

summary(pair_cc_ps$mw)

### create weighted data
pair_cc_ps_svy <- svydesign(ids = ~1,
                                  data = pair_cc_ps,
                                  weights = ~mw)

### weighted table one
table_one_weighted <- svyCreateTableOne(vars = cov_vector3,
                                        strata = "exposure1",
                                        data = pair_cc_ps_svy,
                                        test = FALSE)

table_one_weighted_sav <- print(table_one_weighted, showAllLevels = TRUE, smd = TRUE,
                                nonnormal = nonnorm_vec,
                                formatOptions = list(big.mark = ","))

write.csv(table_one_weighted_sav, "./output/weighted_descriptives/table_one_weighted.csv")

#### Propensity score overlap weight
## Overlap weight
pair_cc_ps$ow <- (pair_cc_ps$ps_assign * (1 - pair_cc_ps$ps_assign)) / pair_cc_ps$ps_assign

summary(pair_cc_ps$ow)
summary(pair_cc_ps$ow[pair_cc_ps$exposure1=="unexposed"])
summary(pair_cc_ps$ow[pair_cc_ps$exposure1!="unexposed"])

## Weighted data
pair_cc_psSvyOw <- svydesign(ids = ~ 1, data = pair_cc_ps, weights = ~ ow)

## Construct a table (This is a bit slow.)
tabWeightedOw <- svyCreateTableOne(vars = cov_vector3, strata = "exposure1", data = pair_cc_psSvyOw, test = FALSE)
## Show table with SMD
print(tabWeightedOw, smd = TRUE)

tabWeightedOw_sav <- print(tabWeightedOw, showAllLevels = TRUE, smd = TRUE,
                                       nonnormal = nonnorm_vec,
                                       formatOptions = list(big.mark = ","))

write.csv(tabWeightedOw_sav, "./output/weighted_descriptives/tabWeightedOw_sav.csv")

################################################################################
#####              inverse probability of treatment weighting              ##### 
################################################################################

#### manual IPW calculation
# Use the 1/PS for the exposed and 1/(1-PS) for the unexposed to create a pseudo-population

pair_cc_manual_ipw <- pair_cc_ps %>% 
  mutate(exposure1 = as.numeric(exposure1)-1) %>% 
  mutate(ipw0 = ifelse(exposure1==0,
                       1/(1 - ps_exp1),0)) %>% 
  mutate(ipw1 = ifelse(exposure1==1,
                       1/(ps_exp1),0)) %>% 
  mutate(ipw = ifelse(ipw0!=0,ipw0,ipw1)) #%>% dplyr::select(ps_exp1, ipw0, ipw1, ipw) #%>% 
#  mutate(ipw_log = log(ipw))

summary(pair_cc_manual_ipw$ipw0)
summary(pair_cc_manual_ipw$ipw1)
summary(pair_cc_manual_ipw$ipw)

ipwplot(pair_cc_manual_ipw$ipw)


summary(pair_cc_manual_ipw$ipw)
#summary(pair_cc_manual_ipw$ipw_log)

### created IPTW dataframe
pair_cc_manual_ipw_svy <- svydesign(ids = ~1,
                             data = pair_cc_manual_ipw,
                             weights = ~ipw)

#pair_cc_manual_ipw_log_svy <- svydesign(ids = ~1,
#                                    data = pair_cc_manual_ipw,
#                                    weights = ~ipw_log)

### create IPTW table one
pair_cc_manual_ipw_tabone <- svyCreateTableOne(vars = cov_vector3,
                                               strata = "exposure1",
                                               data = pair_cc_manual_ipw_svy,
                                               test = FALSE)

pair_cc_manual_ipw_tabone_sav <- print(pair_cc_manual_ipw_tabone, showAllLevels = TRUE, smd = TRUE,
                                       nonnormal = nonnorm_vec,
                                       formatOptions = list(big.mark = ","))

#pair_cc_manual_ipw_log_tabone <- svyCreateTableOne(vars = cov_vector,
#                                                   strata = "exposure1",
#                                                   data = pair_cc_manual_ipw_log_svy,
#                                                   test = FALSE)
#
#pair_cc_manual_ipw_log_tabone_sav <- print(pair_cc_manual_ipw_log_tabone, showAllLevels = TRUE, smd = TRUE,
#                                           nonnormal = nonnorm_vec,
#                                           formatOptions = list(big.mark = ","))

write.csv(pair_cc_manual_ipw_tabone_sav, "./output/temp_output/pair_cc_manual_ipw_tabone_sav.csv")

#write.csv(pair_cc_manual_ipw_log_tabone_sav, "./output/temp_output/pair_cc_manual_ipw_log_tabone_sav.csv")


#### IPW package
pair_cc_ipw <- pair_cc_ps %>% 
  mutate(exposure1 = as.numeric(exposure1)) %>% #,
#         sex_dv_t0 = as.numeric(sex_dv_t0),
#         non_white_t0 = as.numeric(non_white_t0)) %>% 
  mutate(exposure1 = exposure1-1)#,
#         male_t0 = sex_dv_t0-1)#,
#         non_white_t0 = ifelse(non_white_t0==3,0,1)) #%>% 
#  drop_na(everything())
##  dplyr::select(male_t0) %>% 
##  group_by(male_t0) %>% 
##  summarise(n=n())
#
weight <- ipwpoint(exposure = exposure1, family = "binomial", link = "logit",
#                   numerator =~ 1,
                   denominator = ~         sex_dv_t0 +
  age_dv_t0 +
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
                   trunc = 0.01, 
                   data = as.data.frame(pair_cc_ipw))
weight
pair_cc_ipw$.ipw0 <-  weight$weights.trunc
summary(pair_cc_ipw$.ipw0)

pair_cc_ipw$.ipw0_log <-  log(pair_cc_ipw$.ipw0)
summary(pair_cc_ipw$.ipw0_log)


### plot weights
ipw_plot_cc <- ipwplot(pair_cc_ipw$.ipw0)#, timevar = NULL, binwidth = NULL, logscale = TRUE,
#        xlab = NULL, ylab = NULL, main = "", ref = TRUE, ...)

ipwplot(pair_cc_ipw$.ipw0_log)

### created IPTW dataframe
pair_cc_ipw_svy <- svydesign(ids = ~1,
                            data = pair_cc_ipw,
                            weights = ~.ipw0)

pair_cc_ipw_log_svy <- svydesign(ids = ~1,
                             data = pair_cc_ipw,
                             weights = ~.ipw0_log)
### create IPTW table one
table_one_ipw <- svyCreateTableOne(vars = cov_vector,
                                        strata = "exposure1",
                                        data = pair_cc_ipw_svy,
                                        test = FALSE)

table_one_ipw_sav <- print(table_one_ipw, showAllLevels = TRUE, smd = TRUE,
                                nonnormal = nonnorm_vec,
                                formatOptions = list(big.mark = ","))

table_one_ipw_log <- svyCreateTableOne(vars = cov_vector,
                                   strata = "exposure1",
                                   data = pair_cc_ipw_log_svy,
                                   test = FALSE)

table_one_ipw_log_sav <- print(table_one_ipw_log, showAllLevels = TRUE, smd = TRUE,
                           nonnormal = nonnorm_vec,
                           formatOptions = list(big.mark = ","))

write.csv(table_one_ipw_sav, "./output/ipw_descriptives/table_one_ipw.csv")

write.csv(table_one_ipw_log_sav, "./output/temp_output/table_one_ipw_log_sav.csv")

### WeightIt IPTW

weight_weightit <- WeightIt::weightit(exposure1 ~ sex_dv_t0 +
                                        age_dv_t0 +
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
                                        sf12pcs_dv_t0 ,
                                     data = pair_cc_ipw, 
                                     stabilize = TRUE,
                                     estimand = "ATE",  # Find the ATE
                                     method = "ps")  # Build weights with propensity scores
weight_weightit

summary(weight_weightit)
summary(weight_weightit$weights)
density(weight_weightit$weights)

pair_cc_weightit <- pair_cc_ipw
pair_cc_weightit$weightit_ipw <- weight_weightit$weights

pair_cc_weightit %>% group_by(exposure1) %>% summarise(n=sum(weightit_ipw))

### created IPTW dataframe
pair_cc_weightit_ipw_svy <- svydesign(ids = ~1,
                             data = pair_cc_weightit,
                             weights = ~weightit_ipw)

### create IPTW table one
table_one_weightit_ipw <- svyCreateTableOne(vars = cov_vector3,
                                   strata = "exposure1",
                                   data = pair_cc_weightit_ipw_svy,
                                   test = FALSE)

table_one_weightit_ipw_sav <- print(table_one_weightit_ipw, showAllLevels = TRUE, smd = TRUE,
                           nonnormal = nonnorm_vec,
                           formatOptions = list(big.mark = ","))

write.csv(table_one_weightit_ipw_sav, "./output/temp_output/table_one_weightit_ipw_sav.csv")

#### write working files -------------------------------------------------------
write_rds(df_matched, "working_data/pair_cc_matched.rds") # PS matched data
write_rds(pair_cc_ps, "working_data/pair_cc_ps.rds") # analytic df with PS
#write_rds(pair_cc_ipw, "working_data/pair_cc_ipw.rds") #
write_rds(matchit_df, "working_data/matchit_df.rds") # analytic df with IPTW from MatchIt package
write_rds(pair_cc_weightit, "working_data/pair_cc_weightit.rds") # analytic df with IPTW from WeightIt package


################################################################################
#####                  assess balance pre/post-matching                   ######
################################################################################

#### create df for matched smd values ------------------------------------------
table_one_matched_smd <- data.frame(ExtractSmd(table_one_matched))
table_one_matched_smd <- table_one_matched_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "PS matched")

#### create df for PS weighted smd values -----------------------------------------
table_one_weighted_smd <- data.frame(ExtractSmd(table_one_weighted))
table_one_weighted_smd <- table_one_weighted_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "PS weighted")

#### create df for IPW weighted smd values -----------------------------------------
#table_one_ipw_smd <- data.frame(ExtractSmd(table_one_ipw))
#table_one_ipw_smd <- table_one_ipw_smd %>% 
#  rownames_to_column("var") %>% # Apply rownames_to_column
#  rename("smd" = "X1.vs.2") %>% 
#  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
#         matched = "ipw")

#### create df for IPW weighted smd values -----------------------------------------
### MatchIt --------
table_one_MatchIt_smd <- data.frame(ExtractSmd(table_one_matchit))
table_one_MatchIt_smd <- table_one_MatchIt_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "IPTW (MatchIt)")

### WeightIt --------
#table_one_WeightIt_smd <- data.frame(ExtractSmd(table_one_weightit_ipw))
#table_one_WeightIt_smd <- table_one_WeightIt_smd %>% 
#  rownames_to_column("var") %>% # Apply rownames_to_column
#  rename("smd" = "X1.vs.2") %>% 
#  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
#         matched = "IPTW (WeightIt)")


#### bind to unmatched smd df --------------------------------------------------
assess_matching_balance <- table_one_unmatched_smd %>% 
  bind_rows(table_one_matched_smd, 
            table_one_MatchIt_smd) %>% 
  mutate(smd_flag = ifelse(smd>0.2,1,0)) %>% 
  rename("method" = "matched")

### plot
unemp_t1_balance_plot <- assess_matching_balance %>% 
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

tiff("./output/unemp_t1_balance_plot.tiff", width = 960)
unemp_t1_balance_plot
dev.off()

#### number of individuals -----------------------------------------------------
iptw_df %>% dplyr::select(pidp) %>% unique() %>% nrow()

