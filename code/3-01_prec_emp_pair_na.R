################################################################################

# Precarious employment and health - Understanding Society
# 3-01 - Missing data 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Checks item missingness across survey waves 


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
  filter(jr_analytic_var==1) %>% 
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

