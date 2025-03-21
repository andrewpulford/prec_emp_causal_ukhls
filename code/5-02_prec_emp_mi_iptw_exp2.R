################################################################################

# Precarious employment and health - Understanding Society
# 5-02 - IPTW for multiple imputed outcome models  - alternate exposure definition
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a)  calculates IPTW for final MI outcome models
# (b)  produces love plot to show balance between treatment groups
# (c)  produces Table 1 descriptives


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
imputed_data <- readRDS("./working_data/mi/imputed_data_exp2.rds")

## check forNAs in imputed data
sapply(complete(imputed_data,"long"), function(x) sum(is.na(x)))
sapply(complete(imputed_data,"long"), function(x) sum(x=="missing"))


#### prepare data -------------------------------------------------------------- 


################################################################################
#####               inverse probability of treatment weighting             #####
################################################################################

start_time <- Sys.time()
weightit_df <- weightthem(exp2_bin ~
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

### save so you don't have to re-run weighting
write_rds(weightit_df,"./working_data/mi/weightit_df_exp2.rds")

