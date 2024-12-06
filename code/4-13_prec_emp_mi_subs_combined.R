################################################################################

# Precarious employment and health - Understanding Society
# 4-13 - IPTW multiple imputed analytic sample - sub-group analysis final dataset
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a)  produces final dataframe for sub-group analyses
# (b)  


################################################################################

## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
#library(mice) # for multiple imputation
#library(ggmice) # for plotting MI
#library(glmmTMB) # for multi-level modelling
#library(broom.mixed) # for tidying glmmTMB models into df's
#library(cobalt) # Covariate Balance Tables and Plots
#library(MatchThem) # to perform propensity score weighting within each imputation
#library(survey)
#library(glmmTMB) # for multi-level modelling (faster than lme4)
#remotes::install_github("ngreifer/MatchIt")
#library(miceadds) # for working with imputed dataset

################################################################################
#####                         load and prepare data                        #####
################################################################################

female <- read_csv("./output/mi/weighted_outcomes/mi_dr_iptw_df_sex_f.csv") %>% mutate(group="sex")
male <- read_csv("./output/mi/weighted_outcomes/mi_dr_iptw_df_sex_m.csv") %>% mutate(group="sex")
younger <- read_csv("./output/mi/weighted_outcomes/mi_dr_iptw_df_age_younger.csv") %>% mutate(group="age")
older <- read_csv("./output/mi/weighted_outcomes/mi_dr_iptw_df_age_older.csv") %>% mutate(group="age")
rel_pov <- read_csv("./output/mi/weighted_outcomes/mi_dr_iptw_df_rel_pov_y.csv") %>% mutate(group="relative poverty")
rel_pov_no <- read_csv("./output/mi/weighted_outcomes/mi_dr_iptw_df_rel_pov_n.csv") %>% mutate(group="relative poverty")

combined_df <- female %>% 
  bind_rows(male, younger, older, rel_pov, rel_pov_no) %>% 
  dplyr::select(-c(...1,)) %>% 
  dplyr::select(group, sub_group, everything()) %>% 
  arrange(group, est_type ,outcome) %>% 
  mutate(estimate = ifelse(est_type=="OR", exp(estimate),estimate),
         lci = ifelse(est_type=="OR", exp(lci),lci),
         uci= ifelse(est_type=="OR", exp(uci),uci))

