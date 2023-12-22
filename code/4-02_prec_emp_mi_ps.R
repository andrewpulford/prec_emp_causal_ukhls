################################################################################

# Precarious employment and health - Understanding Society
# 4-02 - Propensity score matched multiple imputed analytic sample 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Checks item missingness across survey waves 
# (b) create final MI df 


################################################################################

## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
library(mice) # for multiple imputation
library(ggmice) # for plotting MI
library(glmmTMB) # for multi-level modelling
library(broom.mixed) # for tidying glmmTMB models into df's
library(cobalt) # Covariate Balance Tables and Plots
library(MatchThem) # to perform propensity score weighting within each imputation
library(survey)
library(glmmTMB) # for multi-level modelling (faster than lme4)
#remotes::install_github("ngreifer/MatchIt")

################################################################################
#####                         load and prepare data                        #####
################################################################################


#### read in variable vectors --------------------------------------------------
source("./look_ups/variable_vectors.r")

##### load original data without imoutation ------------------------------------
#mi_subset2 <-  readRDS("./working_data/mi/mi_subset2.rds")

#### load imputed data --------------------------------------------------------
imputed_data <- readRDS("./working_data/mi/imputed_data.rds")

#### prepare data -------------------------------------------------------------- 


################################################################################
#####                         propensity score matching                    #####
################################################################################

start_time <- Sys.time()
matchit_df <- matchthem(exp1_bin ~
                         sex_dv_t0 +
                         age_dv_t0 +
                         non_white_t0 +
                         marital_status_t0 +
                         dep_child_bin_t0 +
                         degree_bin_t0 +
                         gor_dv_t0 +
                         sic2007_section_lab_t0 +
                         soc2000_major_group_title_t0 +
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
                       data = imputed_data,
                       method = "nearest", 
                       ratio=3,
                       distance = "glm",
                       estimand = "ATT")
end_time <- Sys.time()
end_time - start_time

matchit_df
summary(matchit_df)

#with(matchit_df, summary(as.data.frame(mget(ls()))))

#imputed_data$weights_ps <- matchit_df$weights

#with(imputed_data, sum(imputed_data$weights_ps))

test <- bal.tab(matchit_df, un = TRUE, 
                binary = "std", continuous = "std")
test2 <- test$Balance.Across.Imputations

## probably don't need these....
bal.plot(matchit_df, which.imp = 1, 
         var.name = "sex_dv_t0", 
         which = "both")
bal.plot(matchit_df, which.imp = 1, 
         var.name = "age_dv_t0", 
         which = "both")
bal.plot(matchit_df, which.imp = 1, 
         var.name = "non_white_t0", 
         which = "both")

## create love plot to visualise balance between unmatched and matched data across MIs
love.plot(matchit_df, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unmatched","Matched"))
# var.names() - to clean up names

################################################################################
####  #           propensity score matched double-robust model             #####
################################################################################

#### SF-12 PCS -----------------------------------------------------------------

matchit_mods_pcs <- with(matchit_df,
                     glmmTMB(sf12pcs_dv_t1 ~
                               exposure1 +
                               sf12pcs_dv_t0 +
                               sex_dv_t0 +
                               age_dv_t0 +
                               age_dv_t1 +
                               non_white_t0 +
                               marital_status_t0 +
                               dep_child_bin_t0 +
                               degree_bin_t0 +
                               gor_dv_t0 +
                               sic2007_section_lab_t0 +
                               soc2000_major_group_title_t0 +
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
                               (1|pidp)))

matchit_pooled_pcs <- pool(matchit_mods_pcs)

matchit_pooled_pcs_df <- data.frame(summary(matchit_pooled_pcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)


#### SF-12 MCS -----------------------------------------------------------------

matchit_mods_mcs <- with(matchit_df,
                         glmmTMB(sf12mcs_dv_t1 ~
                                   exposure1 +
                                   sf12mcs_dv_t0 +
                                   sex_dv_t0 +
                                   age_dv_t0 +
                                   age_dv_t1 +
                                   non_white_t0 +
                                   marital_status_t0 +
                                   dep_child_bin_t0 +
                                   degree_bin_t0 +
                                   gor_dv_t0 +
                                   sic2007_section_lab_t0 +
                                   soc2000_major_group_title_t0 +
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
                                   (1|pidp)))

matchit_pooled_mcs <- pool(matchit_mods_mcs)

matchit_pooled_mcs_df <- data.frame(summary(matchit_pooled_mcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)

#### Self-rated health ---------------------------------------------------------

matchit_mods_srh <- with(matchit_df,
                         glmmTMB(srh_bin_t1 ~
                                   exposure1 +
                                   srh_bin_t0 +
                                   sex_dv_t0 +
                                   age_dv_t0 +
                                   age_dv_t1 +
                                   non_white_t0 +
                                   marital_status_t0 +
                                   dep_child_bin_t0 +
                                   degree_bin_t0 +
                                   gor_dv_t0 +
                                   sic2007_section_lab_t0 +
                                   soc2000_major_group_title_t0 +
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
                                   (1|pidp),family = binomial(link="logit")))

matchit_pooled_srh <- pool(matchit_mods_srh)

matchit_pooled_srh_df <- data.frame(summary(matchit_pooled_srh, conf.int = TRUE))

matchit_pooled_srh_df <- matchit_pooled_srh_df %>% 
  mutate(estimate = exp(estimate)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

#### GHQ-12 caseness -----------------------------------------------------------

matchit_mods_ghq <- with(matchit_df,
                         glmmTMB(ghq_case4_t1 ~
                                   exposure1 +
                                   ghq_case4_t0 +
                                   sex_dv_t0 +
                                   age_dv_t0 +
                                   age_dv_t1 +
                                   non_white_t0 +
                                   marital_status_t0 +
                                   dep_child_bin_t0 +
                                   degree_bin_t0 +
                                   gor_dv_t0 +
                                   sic2007_section_lab_t0 +
                                   soc2000_major_group_title_t0 +
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
                                   (1|pidp),family = binomial(link="logit")))

matchit_pooled_ghq <- pool(matchit_mods_ghq)

matchit_pooled_ghq_df <- data.frame(summary(matchit_pooled_ghq, conf.int = TRUE))

matchit_pooled_ghq_df <- matchit_pooled_ghq_df %>% 
  mutate(estimate = exp(estimate)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

