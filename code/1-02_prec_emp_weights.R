################################################################################

# Precarious employment and health - Understanding Society
# 1-02 - Create weight spine for causal analysis
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Creates weight spines for analytic samples a and b

################################################################################

## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(readxl) # for reading excel files
library(foreign) # for reading SPSS files
library(tidyverse) # all kinds of stuff 
library(janitor) # cleaning up


################################################################################
#####                         load and prepare data                        #####
################################################################################

### load in master raw file
master_raw1 <- readRDS("./raw_data/master_raw1.rds")



### select pipd, wv_n, weights, psu, strata
weight_spine <- master_raw1 %>% 
  dplyr::select(pidp, wv_n, indinus_lw, indinub_lw, indinui_lw, strata, psu) %>% 
  filter(wv_n!=1) %>% # take out wave 1 as not needed
  rename("wv_end" = "wv_n") %>% 
  # convert into long format
  pivot_longer(cols = c(3:5), names_to = "wt_name", values_to = "wt_value") %>%
  # create flag for weight use with paired waves
  mutate(wt_flag_pair = ifelse(wv_end==2 & wt_name == "indinus_lw",1,
                   ifelse(wv_end%in%c(3:6) & wt_name == "indinub_lw",1,
                   ifelse(wv_end%in%c(7:10) & wt_name == "indinui_lw",1,0)))) %>% 
  # create flag for weight use with three-wave data
  mutate(wt_flag_trio = ifelse(wv_end==3 & wt_name == "indinus_lw",1,
                               ifelse(wv_end%in%c(4:7) & wt_name == "indinub_lw",1,
                                      ifelse(wv_end%in%c(8:10) & wt_name == "indinui_lw",1,0)))) %>% 
  mutate(wt_valid_flag = ifelse(wt_value!=0,1,0))

### create paired wave weight spine
weight_spine_pair <- weight_spine %>%   filter(wt_flag_pair==1) %>% 
  rename("wv_n_t1"="wv_end")

table(weight_spine_pair$wt_flag_pair)

write_rds(weight_spine_pair, "./look_ups/weight_spine_pair.rds")

### create three-wave weight spine
weight_spine_trio <- weight_spine %>%   filter(wt_flag_trio==1) %>% 
  rename("wv_n_t2"="wv_end")

table(weight_spine_trio$wt_flag_trio)

write_rds(weight_spine_trio, "./look_ups/weight_spine_trio.rds")
