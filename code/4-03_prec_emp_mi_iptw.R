################################################################################

# Precarious employment and health - Understanding Society
# 4-03 - IPTW multiple imputed analytic sample 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a)  
# (b)  


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

## check forNAs in imputed data
sapply(complete(imputed_data,"long"), function(x) sum(is.na(x)))


#### prepare data -------------------------------------------------------------- 


################################################################################
#####               inverse probability of treatment weighting             #####
################################################################################

start_time <- Sys.time()
weightit_df <- weightthem(exp1_bin ~
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
                        stabilize = TRUE,
                        estimand = "ATE",  
                        method = "ps")
end_time <- Sys.time()
end_time - start_time

weightit_df
summary(weightit_df)

#with(weightit_df, summary(as.data.frame(mget(ls()))))

#imputed_data$weights_ps <- weightit_df$weights

#with(imputed_data, sum(imputed_data$weights_ps))

test <- bal.tab(weightit_df, un = TRUE, 
                binary = "std", continuous = "std")
test2 <- test$Balance.Across.Imputations

## probably don't need these....
bal.plot(weightit_df, which.imp = 1, 
         var.name = "sex_dv_t0", 
         which = "both")
bal.plot(weightit_df, which.imp = 1, 
         var.name = "age_dv_t0", 
         which = "both")
bal.plot(weightit_df, which.imp = 1, 
         var.name = "non_white_t0", 
         which = "both")

## create love plot to visualise balance between unmatched and matched data across MIs
love.plot(weightit_df, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
# var.names() - to clean up names

################################################################################
#####                               Descriptives                           #####
################################################################################

### leave for now
# use createtableone and leave SMDs for next script
#table_one_ps <- with(weightit_df, CreateTableOne(
#  vars = cov_vector,
#  data=as.data.frame(mget(ls())), 
#  factorVars=catVars_short_vec,
#  strata = "exposure1", test =TRUE))


################################################################################
####  #           propensity score matched double-robust model             #####
################################################################################

#### SF-12 PCS -----------------------------------------------------------------

weightit_mods_pcs <- with(weightit_df,
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

weightit_pooled_pcs <- pool(weightit_mods_pcs)

weightit_pooled_pcs_df <- data.frame(summary(weightit_pooled_pcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)


#### SF-12 MCS -----------------------------------------------------------------

weightit_mods_mcs <- with(weightit_df,
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

weightit_pooled_mcs <- pool(weightit_mods_mcs)

weightit_pooled_mcs_df <- data.frame(summary(weightit_pooled_mcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)

#### Self-rated health ---------------------------------------------------------

weightit_mods_srh <- with(weightit_df,
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

weightit_pooled_srh <- pool(weightit_mods_srh)

weightit_pooled_srh_df <- data.frame(summary(weightit_pooled_srh, conf.int = TRUE))

weightit_pooled_srh_df <- weightit_pooled_srh_df %>% 
  mutate(estimate = exp(estimate)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

#### GHQ-12 caseness -----------------------------------------------------------

weightit_mods_ghq <- with(weightit_df,
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

weightit_pooled_ghq <- pool(weightit_mods_ghq)

weightit_pooled_ghq_df <- data.frame(summary(weightit_pooled_ghq, conf.int = TRUE))

weightit_pooled_ghq_df <- weightit_pooled_ghq_df %>% 
  mutate(estimate = exp(estimate)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

################################################################################
#####                   combine estimates into single df                   #####
################################################################################

pcs_df <- weightit_pooled_pcs_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "SF-12 PCS",
         est_type = "coefficient")

mcs_df <- weightit_pooled_mcs_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "SF-12 MCS",
         est_type = "coefficient")

srh_df <- weightit_pooled_srh_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "Poor self-rated health",
         est_type = "OR")

ghq_df <- weightit_pooled_ghq_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "GHQ-12 caseness (4+)",
         est_type = "OR")

combined_df <- pcs_df %>% 
  bind_rows(mcs_df, srh_df, ghq_df) %>% 
  dplyr::select(c(outcome,	est_type,	estimate,	std.error,	p.value,	lci,	uci))

write.csv(combined_df, "./output/mi/weighted_outcomes/mi_dr_iptw_df.csv")


##################


sapply(complete(imputed_data,2), function(x) sum(is.na(x)))
## 2 NA for sex - try to filter these and inspect
