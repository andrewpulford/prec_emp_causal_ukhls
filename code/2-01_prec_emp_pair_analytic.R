################################################################################

# Precarious employment and health - Understanding Society
# 2-01 - Create complete case analytic sample 
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


####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")

## extra vars for creating exposure vars
extra_vars <- c("jbstat_t1", "nunmpsp_dv_t1",
                "sf12mcs_dv_t0", "sf12mcs_dv_t1")

#### load in variable spine and create vector ----------------------------------
#vars <- readRDS("./variables/vars_master.rds")
#  
#analytic_vars <- vars %>% 
#  filter(type %in% c("exp","cov")) %>% 
#  mutate(variable_t0 = paste0(variable,"_t0"),
#         variable_t1 = paste0(variable,"_t1"))
#
#temp1 <- analytic_vars$variable_t0
#temp2 <- analytic_vars$variable_t1
#
#analytic_vars_vec <- c(temp1,temp2)
#rm(temp1,temp2)

####load eligible cases --------------------------------------------------------
pair_no_att <- readRDS("./working_data/pair_no_att.rds") #%>% 
#  dplyr::select(pidp, all_of(c(id_wt_vector, 
#                         cov_vector, cov_vector2, 
#                         outcome_vector2,
#                         extra_vars))) %>% 
#  dplyr::select(-all_of(extra_vars))

################################################################################
#####                             create NAs df                            #####
################################################################################

#missing_vector <- c("missing", "Missing",
#                    "inapplicable", "proxy",
#                    "refusal", 
#                    "Only available for IEMB", 
#                    "Not available for IEMB",
#                    "don't know")
#

pair_no_att_na <- pair_no_att %>% 
#  mutate(across(everything(), as.character)) %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(is.na(.x),1,0)))


################################################################################
#####                       item missing descriptives                      #####
################################################################################

n_row <- nrow(pair_no_att_na)

pair_no_att_na <- pair_no_att_na %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(.x))) %>% 
  pivot_longer(cols=everything(), names_to = "variable", values_to = "n_NA") %>% 
  mutate(pc_NA = n_NA/n_row*100) 


################################################################################
#####                     create final complete case df                    #####
################################################################################

## sort out baseline sf-12 vars to ID missing vars
pair_no_att$sf12pcs_dv_t0 <- as.character(pair_no_att$sf12pcs_dv_t0)
pair_no_att$sf12pcs_dv_t0 <- as.numeric(pair_no_att$sf12pcs_dv_t0)

pair_no_att$sf12mcs_dv_t0 <- as.character(pair_no_att$sf12mcs_dv_t0)
pair_no_att$sf12mcs_dv_t0 <- as.numeric(pair_no_att$sf12mcs_dv_t0)



### remove any incomplete cases for CC analysis
pair_cc_analytic <- pair_no_att %>% 
  na.omit()




# check
sapply(pair_cc_analytic, function(x) sum(is.na(x)))


#### VVV retain for now - may need for MI VVV ####
# Exclude:
# missing employment status at t1
# missing outcomes - SRH, GHQ-12, SF-12 PCS, SF-12 MCS 

# create flag for missing employment status at t1
#pair_cc_analytic <- pair_no_att  %>% 
#  mutate(exposure_na = ifelse(jbstat_t1%in%missing_vector,1,0)) 

# number of cases excluded due to missing exposure at t1
#pair_cc_analytic_na <- pair_cc_analytic %>% 
#  summarise(n = sum(exposure_na)) %>% 
#  mutate(na_lab = "exposure missing")

#pair_cc_analytic <- pair_cc_analytic %>%
# exclude missing exposure cases
#  filter(exposure_na==0) %>% 
# create flag for missing outcome at t1
#  mutate(outcome_na = ifelse(srh_dv_t1%in%missing_vector,1,
#                      ifelse(ghq_case4_t1%in%missing_vector,1,
#                      ifelse(sf12mcs_dv_t1%in%missing_vector,1,
#                      ifelse(sf12mcs_dv_t1%in%missing_vector,1,0))))) 

# number of cases excluded due to missing outcome at t1
#temp <- pair_cc_analytic %>% 
#  summarise(n = sum(outcome_na)) %>% 
#  mutate(na_lab = "outcome missing")

# bind together with exposure na
#pair_cc_analytic_na <- pair_cc_analytic_na %>% 
#  bind_rows(temp)

# exclude missing outcome cases
#pair_cc_analytic <- pair_cc_analytic %>% 
#  filter(outcome_na==0)



################################################################################
#####                               save df's                              #####
################################################################################

## paired eligible complete case df
write_rds(pair_no_att, "./working_data/pair_no_att.rds")


## paired analytic complete case df
write_rds(pair_cc_analytic, "./working_data/cc/pair_cc_analytic.rds")

## exclusions
write_rds(pair_no_att_na, "./working_data/cc/pair_no_att_na.rds")

###### check for GHQ-12 scores

pair_cc_analytic %>% 
  ggplot(aes(x = as.numeric(scghq2_dv_t0))) +
  geom_histogram() +
  facet_wrap(~exposure1)


pair_cc_analytic %>% 
  ggplot(aes(x = as.numeric(scghq2_dv_t1))) +
  geom_histogram() +
  facet_wrap(~exposure1)

table(pair_cc_analytic$scghq2_dv_t0)
table(pair_cc_analytic$scghq2_dv_t1)


################################################################################
##### scrapbook #####
################################################################################

test2 <- pair_no_att %>% filter(pidp=="816041487")

test2_na <- test2 %>% 
  #  mutate(across(everything(), as.character)) %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(is.na(.x),1,0)))



