################################################################################

# Precarious employment and health - Understanding Society
# 5-04 - IPTW multiple imputed outcome model (SF-12 PCS)  - alternate exposure definition
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
weightit_df <- readRDS("./working_data/mi/weightit_df_exp2.rds")

## check forNAs in imputed data
sapply(complete(weightit_df,"long"), function(x) sum(is.na(x)))


#### prepare data -------------------------------------------------------------- 



################################################################################
####  #           propensity score matched double-robust model             #####
################################################################################

#### SF-12 PCS -----------------------------------------------------------------

weightit_mods_pcs <- with(data = weightit_df,
                          exp = glmmTMB(sf12pcs_dv_t1 ~
                                    exposure2 +
                                    sf12pcs_dv_t0 +
                                    sex_bin +
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
                                    rel_pov_bin +
                                    health_t0 +
                                    health_t1 +
                                    # interaction terms
                                    sex_bin*age_dv_t0 +
                                    sex_bin*rel_pov_bin +
                                    age_dv_t0*rel_pov_bin +
                                    (1|pidp)))

weightit_pooled_pcs <- pool(weightit_mods_pcs)

weightit_pooled_pcs_df <- data.frame(summary(weightit_pooled_pcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(term = str_remove(term, "exposure2"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value)))) %>% 
  dplyr::select("outcome", "term",	"est_type",	"estimate",	"std.error",	"p.value",	"lci",	"uci")


write.csv(weightit_pooled_pcs_df, "./working_data/mi/weightit_pooled_pcs_df_exp2.csv")
write_rds(weightit_mods_pcs, "./working_data/mi/weightit_mods_pcs_exp2.rds")
write_rds(weightit_pooled_pcs, "./working_data/mi/weightit_pooled_pcs_exp2.rds")


################## scrapbook #################

#weightit_mods_pcs <- readRDS("./working_data/mi/weightit_mods_pcs_exp2.rds")#
#
#### try avg_comparisons for sub-group analysis
#library("marginaleffects")
#comp.imp <- lapply(weightit_mods_pcs, function(weightit_pooled_pcs) {
#  avg_comparisons(weightit_mods_pcs, 
#                  variables = "exposure2",
#                  by = "sex_bin")
#})
#
#pooled.comp <- mice::pool(comp.imp, dfcom = Inf)
#
#### try setting weights to zero for one group
#
#
#complete_imp <- complete(weightit_df,action = "long", include = TRUE)
#sapply(complete_imp, function(x) sum(is.na(x)))#
#
#head(complete_imp)
#summary(complete_imp$weights)
#
#
#complete_df <- complete_imp %>% 
#  mutate(weights = ifelse(sex_bin=="Female",weights,0))
#
#
#wiminds_df <- as.mids(complete_df)
#