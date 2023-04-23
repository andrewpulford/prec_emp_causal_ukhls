################################################################################

# Precarious employment and health - Understanding Society
# 3-01 - Create analytic sample 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Checks item missingness across survey waves 
# (b) create final complete case df 


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

#### load in variable spine and create vector ----------------------------------
vars <- readRDS("./variables/vars_master.rds")
  
analytic_vars <- vars %>% 
  filter(type %in% c("exp","cov")) %>% 
  mutate(variable_t0 = paste0(variable,"_t0"),
         variable_t1 = paste0(variable,"_t1"))

temp1 <- analytic_vars$variable_t0
temp2 <- analytic_vars$variable_t1

analytic_vars_vec <- c(temp1,temp2)
rm(temp1,temp2)

####load eligible cases --------------------------------------------------------
pair_cc_eligible <- readRDS("./working_data/pair_cc_eligible.rds")# %>% 
#  dplyr::select(analytic_vars_vec)


################################################################################
#####                             create NAs df                            #####
################################################################################

missing_vector <- c("missing", 
                    "inapplicable", "proxy",
                    "refusal", 
                    "Only available for IEMB", 
                    "Not available for IEMB",
                    "don't know")

pair_cc_eligible_na <- pair_cc_eligible %>% 
  mutate(across(everything(), as.character)) %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(is.na(.x),"missing",.x))) %>% 
  mutate(across(.cols=everything(), 
                .fns= ~ifelse(.x%in%missing_vector,1,0)))


################################################################################
#####                       item missing descriptives                      #####
################################################################################

n_row <- nrow(pair_cc_eligible_na)

pair_cc_eligible_na <- pair_cc_eligible_na %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(.x))) %>% 
  pivot_longer(cols=1:120, names_to = "variable", values_to = "n_NA") %>% 
  mutate(pc_NA = n_NA/n_row*100)

################################################################################
#####                     create final complete case df                    #####
################################################################################

# Exclude:
# missing employment status at t1
# missing outcomes - SRH, GHQ-12, SF-12 PCS, SF-12 MCS 

# create flag for missing employment status at t1
pair_cc_analytic <- pair_cc_eligible  %>% 
  mutate(exposure_na = ifelse(jbstat_t1%in%missing_vector,1,0)) 

# number of cases excluded due to missing exposure at t1
pair_cc_analytic_na <- pair_cc_analytic %>% 
  summarise(n = sum(exposure_na)) %>% 
  mutate(na_lab = "exposure missing")

pair_cc_analytic <- pair_cc_analytic %>%
# exclude missing exposure cases
  filter(exposure_na==0) %>% 
# create flag for missing outcome at t1
  mutate(outcome_na = ifelse(srh_dv_t1%in%missing_vector,1,
                      ifelse(ghq_case4_t1%in%missing_vector,1,
                      ifelse(sf12mcs_dv_t1%in%missing_vector,1,
                      ifelse(sf12mcs_dv_t1%in%missing_vector,1,0))))) 

# number of cases excluded due to missing outcome at t1
temp <- pair_cc_analytic %>% 
  summarise(n = sum(outcome_na)) %>% 
  mutate(na_lab = "outcome missing")

# bind together with exposure na
pair_cc_analytic_na <- pair_cc_analytic_na %>% 
  bind_rows(temp)

# exclude missing outcome cases
pair_cc_analytic <- pair_cc_analytic %>% 
  filter(outcome_na==0)

################################################################################
#####                       create exposure variables                      #####
################################################################################

pair_cc_analytic <- pair_cc_analytic %>%
  # unemployed at t1
  mutate(exposure1 = ifelse(jbstat_t1%in%c("unemployed","Unemployed"),
                            "exposed (unemployed at t1)","unexposed")) %>% 
  # job loss between t0 and t1
  mutate(exposure2 = ifelse(jbstat_t1 %in% c("unemployed","unemployed"), 
                            "exposed (job loss between t0 and t1",
                           ifelse(nunmpsp_dv_t1>0,
                                  "exposed (job loss between t0 and t1",
                                  "unexposed"))) 

################################################################################
#####                               save df's                              #####
################################################################################

## paired analytic complete case df
write_rds(pair_cc_analytic, "./working_data/pair_cc_analytic.rds")

## exclusions
write_rds(pair_cc_analytic_na, "./working_data/pair_cc_analytic_na.rds")

