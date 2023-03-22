################################################################################

# Persistent precarious employment and health - Understanding Society
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

`%notin%` <- Negate(`%in%`)


################################################################################
#####                         load and prepare data                        #####
################################################################################

pair_cc_raw <- readRDS("./working_data/pair_cc_raw.rds")


################################################################################
#####                             job retention                            #####
################################################################################


#### apply inclusion criteria ---------------------------------------------------
# in employment t0
# valid employment response at t1
# 16-64 t0
# not in full-time education at t0 or t1
# not retired at t0 or t1


### create df and exclusions df
## in employment t0
# include
in_emp_t0 <- pair_cc_raw %>% 
  filter(employ_t0=="yes") 

# exclude
temp_df <- pair_cc_raw %>% filter(employ_t0!="yes") 
pair_cc_exc <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "not in employment t0")

## valid employment response at t1
# include
valid_emp_t1 <- in_emp_t0 %>% 
  filter(employ_t1=="yes" | employ_t1=="no") 

# exclude
temp_df <- in_emp_t0 %>%   filter(employ_t1%notin%c("yes", "no"))
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "no valid employment response at t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)


## 16-64 t0
# include
working_age <- valid_emp_t1 %>% 
  filter(age_dv_t0>=16 &age_dv_t0<=64)

# exclude
temp_df <- valid_emp_t1 %>% filter(age_dv_t0<16 | age_dv_t0>64) 
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "not working age (16-64 years) at t0")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)


## not in full-time education at t0 or t1
# inclue
not_fe <- working_age %>% 
  filter(fenow_t0!="At college/university" |
           fenow_t1!="At college/university")

# exclude
temp_df <- working_age %>%  filter(fenow_t0=="At college/university" |
                                 fenow_t1=="At college/university")
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "at college/university at t0 or t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)

## not retired at t0 or t1
# include
not_retired <- not_fe %>% 
  filter(retchk_flag_t0!=1 | retchk_flag_t1!=1)

# exclude
temp_df <- not_fe %>% filter(retchk_flag_t0==1 | retchk_flag_t1==1)
temp2 <- temp_df %>% summarise(n=n()) %>% 
  mutate(exc_reason = "retired at t0 or t1")
pair_cc_exc <- pair_cc_exc %>% 
  bind_rows(temp2)


#### final df's
pair_cc_eligible <- not_retired

##save

write_rds(pair_cc_eligible, "./working_data/pair_cc_eligible.rds")
