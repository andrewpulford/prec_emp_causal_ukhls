################################################################################

# Precarious employment and health - Understanding Society
# 1-05 - Create paired analytic sample dataframes
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) creates eligible person-trials dataframes
# (b) creates summary df for excluded cases

################################################################################

## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
library(janitor) # cleaning up
library(naniar) # for missing values

`%notin%` <- Negate(`%in%`)


################################################################################
#####                         load and prepare data                        #####
################################################################################

####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")


## extra vars for 
extra_vars <- c("jbstat_t0", "jbstat_t1","nunmpsp_dv_t0", 
                "nunmpsp_dv_t1", "retchk_flag_t1",
                "scghq2_dv_t0", "scghq2_dv_t1")

#### read in data
pair_cc_raw <- readRDS("./working_data/pair_cc_raw.rds")

weight_spine_pair <- readRDS("./look_ups/weight_spine_pair.rds")

pair_cc_raw <- pair_cc_raw %>% 
  left_join(weight_spine_pair) %>% 
  dplyr::select(all_of(c(id_wt_vector, 
                         cov_vector, cov_vector2, 
                         outcome_vector2,
                         extra_vars))) %>% 
  dplyr::select(-c(psu, strata, wt_name, wt_value))

  

################################################################################
#####                             job retention                            #####
################################################################################


#### apply inclusion criteria ---------------------------------------------------
# 16-64 during study period
# in employment t0
# valid employment response at t1
# not self-employed or undertaking unpaid work for family business at T1 
# not on a government training scheme or apprenticeship at T1
# not in full-time education at t0 or t1
# not retired at t0 or t1
# not on maternity leave at T1
# not caring for family or home at T1
# not long-term sick or disabled at T1
# not other unspecified employment status at t1
# not inapplicable or proxy outcome measure at t0 or t1

### create df and exclusions df
## 16-64 t0
# include
working_age <- pair_cc_raw  %>% 
  filter(age_dv_t0>=16 &age_dv_t1<=64)

# exclude
temp_df <- pair_cc_raw  %>% filter(age_dv_t0<16 | age_dv_t1>64 |
                                     is.na(age_dv_t0) | is.na(age_dv_t1)) 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "not working age (16-64 years)")
pair_cc_exc <- temp2


## in employment t0
# include
in_emp_t0 <-  working_age %>% 
  filter(jbstat_t0=="Paid employment(ft/pt)") 

# exclude
temp_df <- working_age %>% filter(jbstat_t0!="Paid employment(ft/pt)") 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "not in employment t0")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)


################################################################################
#####                       create exposure variables                      #####
################################################################################

pair_eligible <- in_emp_t0 %>%
  # prevention of unemployment at t1
  mutate(exposure1 = ifelse(jbstat_t1%in%c("unemployed","Unemployed"),
                            "unexposed",
                            ifelse(jbstat_t1=="Paid employment(ft/pt)",
                                   "exposed (employed at t1)" ,
                                   "CHECK"))) %>% 
  # prevention of job loss between t0 and t1
  mutate(exposure2 = ifelse(jbstat_t1 %in% c("unemployed","unemployed") | nunmpsp_dv_t1>0, 
                            "unexposed",
                            ifelse(jbstat_t1=="Paid employment(ft/pt)" & nunmpsp_dv_t1==0,
  "exposed (no job loss between t0 and t1", "CHECK")))


################################################################################
#####                               attrition                              #####
################################################################################

## not self-employed or undertaking unpaid work for family business at T1
# include
not_self_emp <- pair_eligible %>% 
  filter(jbstat_t1%notin%c("Self employed", "self employed", 
                           "Unpaid, family business"))

# exclude
temp_df <- pair_eligible %>%  filter(jbstat_t1%in%c("Self employed", 
                                                "self employed", 
                                                "Unpaid, family business"))
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "self-employed or undertaking unpaid work for family business at t1")
pair_cc_att <- temp2


## not government training scheme or apprenticeship at T1 
# include
not_training <- not_self_emp %>% 
  filter(jbstat_t1%notin%c("Govt training scheme", "On apprenticeship"))

# exclude
temp_df <- not_self_emp %>%  filter(jbstat_t1%in%c("Govt training scheme", 
                                                "On apprenticeship"))
temp2 <- temp_df %>%   
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "government training scheme or apprenticeship at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

## not in full-time education at t1
# include
not_fe <- not_training %>% 
  filter(jbstat_t1 %notin% c("full-time student", "Full-time student"))

# exclude
temp_df <- not_training %>%  filter(jbstat_t1 %in% c("full-time student", 
                                                     "Full-time student"))
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "full-time student at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

## not retired at t1
# include
not_retired <- not_fe %>% 
  filter(retchk_flag_t1!=1 & 
           jbstat_t1 %notin% c("Retired","retired"))

# exclude
temp_df <- not_fe %>% filter(retchk_flag_t1==1 | 
                               jbstat_t1 %in% c("Retired","retired")) 
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "retired at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)


## not on maternity leave
# include 
not_mat <- not_retired %>% 
  filter(jbstat_t1 %notin% c("on maternity leave","On maternity leave"))

# exclude
temp_df <- not_retired %>% filter(jbstat_t1 %in% c("on maternity leave",
                                                   "On maternity leave")) 
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "on maternity leave at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

## not caring for family or home
# include 
not_caring <- not_mat %>% 
  filter(jbstat_t1!="Family care or home")

# exclude
temp_df <- not_mat %>% filter(jbstat_t1=="Family care or home") 
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "caring for family or home at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

## not long-term sick or disabled at t1
# include 
not_lts <- not_caring %>% 
  filter(jbstat_t1!="LT sick or disabled")

# exclude
temp_df <- not_caring %>% filter(jbstat_t1=="LT sick or disabled") 
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "long-term sick or disabled at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

# not other unspecified employment status at t1
# include 
not_other <- not_lts %>% 
  filter(jbstat_t1%notin%c("doing something else","Doing something else",
                           "Missing", "missing",                                                     
                           "inapplicable",  "proxy",                                                       
                           "refusal", "don't know",                                                  
                           "Only available for IEMB",                                     
                           "Not available for IEMB"))

# exclude
temp_df <- not_lts %>% filter(jbstat_t1%in%c("doing something else","Doing something else",
                                             "Missing", "missing",                                                     
                                             "inapplicable",  "proxy",                                                       
                                             "refusal", "don't know",                                                  
                                             "Only available for IEMB",                                     
                                             "Not available for IEMB")) 
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "other unspecified employment status at t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

## not inapplicable or proxy outcome measure at t0 or t1
# include
valid_outcomes <- not_other %>% 
  filter(sf12pcs_dv_t0%notin%c("inapplicable/proxy","missing")&
           sf12pcs_dv_t1%notin%c("inapplicable/proxy","missing")&
           sf12mcs_dv_t0%notin%c("inapplicable/proxy","missing")&
           sf12pcs_dv_t1%notin%c("inapplicable/proxy","missing")&
           srh_bin_t0%notin%c("inapplicable/proxy","missing")&
           srh_bin_t1%notin%c("inapplicable/proxy","missing")&
           ghq_case4_t0%notin%c("inapplicable/proxy","missing")&
           ghq_case4_t1%notin%c("inapplicable/proxy","missing"))


# exclude
temp_df <- not_other %>% filter(sf12pcs_dv_t0%in%c("inapplicable/proxy","missing")|
                                  sf12pcs_dv_t1%in%c("inapplicable/proxy","missing")|
                                  sf12mcs_dv_t0%in%c("inapplicable/proxy","missing")|
                                  sf12pcs_dv_t1%in%c("inapplicable/proxy","missing")|
                                  srh_bin_t0%in%c("inapplicable/proxy","missing")|
                                  srh_bin_t1%in%c("inapplicable/proxy","missing")|
                                  ghq_case4_t0%in%c("inapplicable/proxy","missing")|
                                  ghq_case4_t1%in%c("inapplicable/proxy","missing")) 
temp2 <- temp_df %>% 
  group_by(exposure1) %>% 
  summarise(n=n()) %>% 
  mutate(exc_reason = "health outcomes not measured at t0 or t1")
pair_cc_att <- pair_cc_att %>% 
  bind_rows(temp2)

#### final df's
pair_no_att <- valid_outcomes


### recode missing categories as NA

pair_no_att <- pair_no_att %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(.x%in%c("Missing",
                                       "missing",                                                     
                                       "inapplicable",                                                
                                       "proxy",                                                       
                                       "refusal",                                                     
                                       "don't know",                                                  
                                       "Only available for IEMB",                                     
                                       "Not available for IEMB"),NA,.x))) 

# check
sapply(pair_no_att, function(x) sum(is.na(x)))

################################################################################
#####                                 save                                 #####
################################################################################

## eligible/no attrition df for CC and MI analysis
write_rds(pair_no_att, "./working_data/pair_no_att.rds")
write_rds(pair_no_att, "./working_data/pair_mi.rds")

## exclusion criteria
write_rds(pair_cc_exc, "./working_data/pair_cc_exc.rds")
write_csv(pair_cc_exc, "./output/cc/pair_cc_exc.csv")

## attrition criteria
write_rds(pair_cc_att, "./working_data/pair_cc_att.rds")
write_csv(pair_cc_att, "./output/cc/pair_cc_att.csv")

################################################################################
#####                            scrapbook                                 #####
################################################################################

table(pair_no_att$exposure1)
table(pair_no_att$exposure2)

exp1_check <- pair_no_att %>% filter(exposure1=="CHECK")

table(exp1_check$jbstat_t1)
table(exp1_check$nunmpsp_dv_t1, exp1_check$exposure1)

