################################################################################

# Precarious employment and health - Understanding Society
# 4-07 - IPTW multiple imputed outcome model (GHQ-12 caseness) 
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
weightit_df <- readRDS("./working_data/mi/weightit_df.rds")

## check forNAs in imputed data
sapply(complete(weightit_df,"long"), function(x) sum(is.na(x)))


#### prepare data -------------------------------------------------------------- 



################################################################################
#####             propensity score matched double-robust model             #####
################################################################################

#### GHQ-12 caseness ----------------------------------------------------

## fit model to the imputations and pool the results:
weightit_mods_ghq <- with(data = weightit_df, 
                          exp = glmmTMB(ghq_case4_t1 ~
                                          # exposure
                                          exposure1 +
                                          # t0 outcome measure
                                          ghq_case4_t0 +
                                          # t0 covariates
                                          sex_dv_t0 +
                                          age_dv_t0 +
                                          non_white_t0 +
                                          marital_status_t0 +
                                          degree_bin_t0 +
                                          gor_dv_t0 +
                                          sic2007_section_lab_t0 +
                                          soc2000_major_group_title_t0 +
                                          jbft_dv_t0 +
                                          small_firm_t0 +
                                          emp_contract_t0 +
                                          broken_emp_t0 +
                                          j2has_dv_t0 +
                                          rel_pov_t0 +
                                          health_t0 +
                                          # t1 covariates
                                          age_dv_t1 +
                                          marital_status_t1 +
                                          health_t1 +
                                          # interaction terms
                                          sex_dv_t0*age_dv_t0 +
                                          sex_dv_t0*rel_pov_t0 +
                                          age_dv_t0*rel_pov_t0 +
                                          (1|pidp),
                                        family=binomial(link="logit")))

weightit_pooled_ghq <- pool(weightit_mods_ghq)

weightit_pooled_ghq_df <- data.frame(summary(weightit_pooled_ghq, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value)))) %>% 
  dplyr::select("outcome", "term",	"est_type",	"estimate",	"std.error",	"p.value",	"lci",	"uci")


write.csv(weightit_pooled_ghq_df, "./working_data/mi/weightit_pooled_ghq_df.csv")

################################################################################
##### scrapbook
################################################################################


## fit model to the imputations and pool the results:
test <- with(data = weightit_df, 
                          exp = glmmTMB(ghq_case4_t1 ~
                                          # exposure
                                          exposure1 +
                                          # t0 outcome measure
#                                          ghq_case4_t0 +
                                          # t0 covariates
 #                                         sex_dv_t0 +
  #                                        age_dv_t0 +
   #                                       non_white_t0 +
    #                                      marital_status_t0 +
     #                                     degree_bin_t0 +
      #                                    gor_dv_t0 +
       #                                   sic2007_section_lab_t0 +
        #                                  soc2000_major_group_title_t0 +
         #                                 jbft_dv_t0 +
          #                                small_firm_t0 +
           #                               emp_contract_t0 +
            #                              broken_emp_t0 +
             #                             j2has_dv_t0 +
              #                            rel_pov_t0 +
               #                           health_t0 +
                                          # t1 covariates
                #                          age_dv_t1 +
                 #                         marital_status_t1 +
                  #                        health_t1 +
                                          # interaction terms
                   #                       sex_dv_t0*age_dv_t0 +
                    #                      sex_dv_t0*rel_pov_t0 +
                     #                     age_dv_t0*rel_pov_t0 +
                                          (1|pidp), 
                                        family=binomial(link="logit")))

test_pooled <- pool(test)

test_pooled_df <- data.frame(summary(test_pooled, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value)))) %>% 
  dplyr::select("outcome", "term",	"est_type",	"estimate",	"std.error",	"p.value",	"lci",	"uci")

