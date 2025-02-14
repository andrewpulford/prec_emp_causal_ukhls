################################################################################

# Precarious employment and health - Understanding Society
# 6-08 - IPTW multiple imputed outcome models (combined df) - full-time only 
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

################################################################################
#####                         load and prepare data                        #####
################################################################################

weightit_pooled_pcs_df <- read.csv("./working_data/mi/weightit_pooled_pcs_df_ft.csv")
weightit_pooled_mcs_df <- read.csv("./working_data/mi/weightit_pooled_mcs_df_ft.csv")
weightit_pooled_srh_df <- read.csv("./working_data/mi/weightit_pooled_srh_df_ft.csv")
weightit_pooled_ghq_df <- read.csv("./working_data/mi/weightit_pooled_ghq_df_ft.csv")

################################################################################
#####               create single summary df for MI outcomes               #####
################################################################################

mi_iptw_df <- weightit_pooled_pcs_df %>% 
  bind_rows(weightit_pooled_mcs_df,
            weightit_pooled_srh_df,  
            weightit_pooled_ghq_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  #  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, est_type, estimate, std.error, p.value, lci, uci) %>% 
  # change binary outcomes to odds ratios
  mutate(est_type = ifelse(outcome %in% c("Poor self-rated health",
                                         "GHQ-12 caseness"),"OR",est_type)) %>% 
  mutate(estimate = ifelse(outcome %in% c("Poor self-rated health",
                                           "GHQ-12 caseness"), exp(estimate),estimate)) %>% 
  mutate(lci = ifelse(outcome %in% c("Poor self-rated health",
                                          "GHQ-12 caseness"), exp(lci),lci)) %>% 
  mutate(uci = ifelse(outcome %in% c("Poor self-rated health",
                                      "GHQ-12 caseness"), exp(uci),uci))

write.csv(mi_iptw_df, "./output/mi/sensitivity_analyses/full-time/weighted_outcomes/mi_iptw_df_ft.csv")


