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
temp_df <- pair_cc_raw  %>% filter(age_dv_t0<16 | age_dv_t1>64) 
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


### valid employment response at t1
# include
#valid_emp_t1 <- in_emp_t0 %>% 
#  filter(employ_t1=="yes" | employ_t1=="no") 
#
## exclude
#temp_df <- in_emp_t0 %>%   filter(employ_t1%notin%c("yes", "no"))
#temp2 <- temp_df %>% summarise(n=n()) %>% 
#  mutate(exc_reason = "no valid employment response at t1")
#pair_cc_exc <- pair_cc_exc %>% 
#  bind_rows(temp2)

## not self-employed or undertaking unpaid work for family business at T1
# include
not_self_emp <- in_emp_t0 %>% 
  filter(jbstat_t1%notin%c("Self employed", "self employed", 
                           "Unpaid, family business"))

# exclude
temp_df <- in_emp_t0 %>%  filter(jbstat_t1%in%c("Self employed", 
                                                      "self employed", 
                                                      "Unpaid, family business"))
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "self-employed or undertaking unpaid work for family business at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)


## not government training scheme or apprenticeship at T1 
# include
not_training <- not_self_emp %>% 
  filter(jbstat_t1%notin%c("Govt training scheme", "On apprenticeship"))

# exclude
temp_df <- in_emp_t0 %>%  filter(jbstat_t1%in%c("Govt training scheme", 
                                                "On apprenticeship"))
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "government training scheme or apprenticeship at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)

## not in full-time education at t1
# include
not_fe <- not_training %>% 
  filter(jbstat_t1 %notin% c("full-time student", "Full-time student"))

# exclude
temp_df <- not_training %>%  filter(jbstat_t1 %in% c("full-time student", 
                                                       "Full-time student"))
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "full-time student at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)

## not retired at t1
# include
not_retired <- not_fe %>% 
  filter(retchk_flag_t1!=1 & 
           jbstat_t1 %notin% c("Retired","retired"))

# exclude
temp_df <- not_fe %>% filter(retchk_flag_t1==1 | 
                             jbstat_t1 %in% c("Retired","retired")) 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "retired at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)


## not on maternity leave
# include 
not_mat <- not_retired %>% 
  filter(jbstat_t1 %notin% c("on maternity leave","On maternity leave"))

# exclude
temp_df <- not_retired %>% filter(jbstat_t1 %in% c("on maternity leave",
                                              "On maternity leave")) 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "on maternity leave at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)

## not caring for family or home
# include 
not_caring <- not_mat %>% 
  filter(jbstat_t1!="Family care or home")

# exclude
temp_df <- not_mat %>% filter(jbstat_t1=="Family care or home") 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "caring for family or home at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)

## not long-term sick or disabled at t1
# include 
not_lts <- not_caring %>% 
  filter(jbstat_t1!="LT sick or disabled")

# exclude
temp_df <- not_caring %>% filter(jbstat_t1=="LT sick or disabled") 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "long-term sick or disabled at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)

# not other unspecified employment status at t1
# include 
not_other <- not_lts %>% 
  filter(jbstat_t1%notin%c("doing something else","Doing something else"))

# exclude
temp_df <- not_lts %>% filter(jbstat_t1%in%c("doing something else",
                                                "Doing something else")) 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "other unspecified employment status at t1")
pair_cc_exc <- pair_cc_exc %>% 
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
temp_df <- not_other %>% filter(sf12pcs_dv_t0=="inapplicable/proxy"|
                                   sf12pcs_dv_t1=="inapplicable/proxy"|
                                   sf12mcs_dv_t0=="inapplicable/proxy"|
                                   sf12pcs_dv_t1=="inapplicable/proxy"|
                                   srh_bin_t0=="inapplicable/proxy"|
                                   srh_bin_t1=="inapplicable/proxy"|
                                   ghq_case4_t0=="inapplicable/proxy"|
                                   ghq_case4_t1=="inapplicable/proxy") 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "health outcomes not measured at t0 or t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)
## valid weight at t1
# include
#valid_wt <- not_retired %>% 
#  filter(wt_valid_flag==1)

# exclude
#temp_df <- not_retired %>% filter(wt_valid_flag!=1)
#temp2 <- temp_df %>% summarise(n=n()) %>% 
#  mutate(exc_reason = "no valid weight")
#pair_cc_exc <- pair_cc_exc %>% 
#  bind_rows(temp2)


#### final df's
pair_eligible <- valid_outcomes

################################################################################
#####                       create exposure variables                      #####
################################################################################

pair_eligible <- pair_eligible %>%
  # prevention of unemployment at t1
  mutate(exposure1 = ifelse(jbstat_t1%in%c("unemployed","Unemployed"),
                            "unexposed","exposed (employed at t1)")) %>% 
  # prevention of job loss between t0 and t1
  mutate(exposure2 = ifelse(jbstat_t1 %in% c("unemployed","unemployed"), 
                            "exposed (no job loss between t0 and t1",
                            ifelse(nunmpsp_dv_t1==0,
                                   "exposed (no job loss between t0 and t1",
                            "unexposed"))) 

### recode missing categories as NA

pair_eligible <- pair_eligible %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 

# check
sapply(pair_eligible, function(x) sum(is.na(x)))


##save

## eligible df for CC and MI analysis
write_rds(pair_eligible, "./working_data/pair_eligible.rds")

write_rds(pair_cc_exc, "./working_data/pair_cc_exc.rds")

