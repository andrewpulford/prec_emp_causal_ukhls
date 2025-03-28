################################################################################

# Precarious employment and health - Understanding Society
# 4-09 - MI EQ5D/QALY analysis   
# analysis for job retention scheme 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) complete case EQ5D analysis for unweighted data
# (b) complete case EQ5D analysis for propensity score matched data
# (c) complete case EQ5D analysis for IPTW data



################################################################################

## remove any existing objects from global environment
rm(list=ls()) 

## turn off scientific notation
options(scipen = 999)

################################################################################
#####                            install packages                          #####
################################################################################

### use tableone v0.12.0 and survey v4.0 so that svy table works
#remotes::install_version("tableone", version = "0.12.0")
#remotes::install_version("survey", version = "4.0")

library(tidyverse) # all kinds of stuff 
library(mice) # for multiple imputation
library(ggmice) # for plotting MI
library(glmmTMB) # for multi-level modelling
library(broom.mixed) # for tidying glmmTMB models into df's
library(cobalt) # Covariate Balance Tables and Plots
library(MatchThem) # to perform propensity score weighting within each imputation
library(survey)
library(glmmTMB) # for multi-level modelling (faster than lme4)
library(tableone)
#library(doParallel)
#library(foreach)
#cluster = makeCluster(detectCores() - 1) # convention to leave 1 core for OS
#registerDoParallel(cluster)

################################################################################
#####                         load and prepare data                        #####
################################################################################


#### load imputed data ---------------------------------------------------------
weightit_df <- readRDS("./working_data/mi/weightit_df.rds")

sapply(complete(weightit_df,"long"), function(x) sum(is.na(x)))


#### prepare data -------------------------------------------------------------- 

### create complete df of all imputations
weightit_df_complete <- complete(weightit_df, action = 'long', include = FALSE) %>% 
  # rename weights as it throws an error otherwise
  rename("ps_weights" = "weights") %>% 
  mutate(dep_child_bin_t0 = ifelse(dep_child_bin_t0==1,"Yes","No"))

sapply(weightit_df_complete, function(x) sum(is.na(x)))


#### parameter estimates for sf-12 - EQ5D mapping ------------------------------
# source: https://journals.sagepub.com/doi/10.1177/0272989X04265477?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
eq5d_franks <- read.csv("./look_ups/eq5d_franks.csv")

## check variable structure
str(eq5d_franks)

## for some reason minus has converted to question mark in negative parameter estimates
## so replace in string then convert to numeric 
eq5d_franks$parameter_est <- str_replace(eq5d_franks$parameter_est,"\\?", "-")
eq5d_franks$parameter_est <- as.numeric(eq5d_franks$parameter_est)

intercept <- eq5d_franks$parameter_est[eq5d_franks$predictor=="intercept"]
PCS <- eq5d_franks$parameter_est[eq5d_franks$predictor=="PCS"]
MCS <- eq5d_franks$parameter_est[eq5d_franks$predictor=="MCS"]
PCSxPCS <- eq5d_franks$parameter_est[eq5d_franks$predictor=="PCSxPCS"]
MCSxMCS <- eq5d_franks$parameter_est[eq5d_franks$predictor=="MCSxMCS"]
PCSxMCS <- eq5d_franks$parameter_est[eq5d_franks$predictor=="PCSxMCS"]

################################################################################
#####                               function                               #####
################################################################################

eq5d_mapper <- function(data, pcs_var, mcs_var){
  data %>% 
    mutate(pcs_para = (pcs_var-49.9)*PCS,
           mcs_para = (mcs_var-51.5)*MCS,
           pcsxpcs_para = ((pcs_var-49.9)^2)*PCSxPCS,
           mcsxmcs_para = ((mcs_var-51.5)^2)*MCSxMCS,
           pcsxmcs = (pcs_var-49.9)*(mcs_var-51.5),
           pcsxmcs_para = pcsxmcs*PCSxMCS) %>% 
    mutate(eq5d_temp = intercept + pcs_para + mcs_para +
             pcsxpcs_para + mcsxmcs_para + pcsxmcs_para)  %>% 
    dplyr::select(-c(pcs_para,mcs_para,pcsxpcs_para,mcsxmcs_para,pcsxmcs,pcsxmcs_para))
}



################################################################################
#####                            unweighted data                           #####
################################################################################

#### t0 EQ-5D ------------------------------------------------------------------
weightit_df_complete <- eq5d_mapper(data=weightit_df_complete, 
                             pcs_var = weightit_df_complete$sf12pcs_dv_t0,
                             mcs_var = weightit_df_complete$sf12mcs_dv_t0) %>% 
  rename(eq5d_t0=eq5d_temp)

summary(weightit_df_complete$eq5d_t0)
hist(weightit_df_complete$eq5d_t0, breaks = 200)


plot(weightit_df_complete$eq5d_t0,weightit_df_complete$sf12pcs_dv_t0)

weightit_df_complete %>% group_by(exposure1) %>% summarise(mean(eq5d_t0))

boxplot(weightit_df_complete$eq5d_t0 ~ weightit_df_complete$exposure1)

#### t1 EQ-5D ------------------------------------------------------------------
weightit_df_complete <- eq5d_mapper(data=weightit_df_complete, 
                             pcs_var = weightit_df_complete$sf12pcs_dv_t1,
                             mcs_var = weightit_df_complete$sf12mcs_dv_t1) %>% 
  rename(eq5d_t1=eq5d_temp)

summary(weightit_df_complete$eq5d_t1)
hist(weightit_df_complete$eq5d_t1, breaks = 200)


plot(weightit_df_complete$eq5d_t1,weightit_df_complete$sf12pcs_dv_t1)

weightit_df_complete %>% group_by(exposure1) %>% summarise(mean(eq5d_t1))

boxplot(weightit_df_complete$eq5d_t1 ~ weightit_df_complete$exposure1)


#### calculate qalys -----------------------------------------------------------
## this is based on approx 12 months treatment
## QALY at single time-point is equal to EQ-5D score at same time-point
## for analysis we calculate the difference between QALY at t1 and t0
qaly_df <- weightit_df_complete %>% mutate(qaly_t0 = eq5d_t0,
                                          qaly_t1 = eq5d_t1,
                                          qaly_diff = qaly_t1-qaly_t0)


### compare qalys for treatment groups (divide by 25 to pool MI data)
qaly_grouped <- qaly_df %>% group_by(exposure1) %>% 
  summarise(qaly_diff = sum(qaly_diff)/25)

### total number of participants by treatment group
n_grouped <- qaly_df %>% group_by(exposure1) %>% 
  summarise(n = n()/25)

## total number in treatment group
n_treated <- sum(qaly_df$exposure1=="exposed (employed at t1)")/25

## calculate average qalys per person per group to get comparable values
qaly_grouped <- qaly_grouped %>% 
  left_join(n_grouped) %>% mutate(qaly_person=qaly_diff/n)

### calculate QALY gain
## QALY gain per treated person
qaly_gain_person <-  qaly_grouped$qaly_person[qaly_grouped$exposure1=="exposed (employed at t1)"] - qaly_grouped$qaly_person[qaly_grouped$exposure1=="unexposed"]

## total for treated for treatment group
qualy_gain_total <- qaly_gain_person*n_treated

### calculate treatment benefit (based on £70k per QALY - UK Govt Green Book)
# https://www.gov.uk/government/publications/the-green-book-appraisal-and-evaluation-in-central-government/the-green-book-2020#valuation-of-costs-and-benefits
benefit <-  qaly_gain_person*70000

### calculate treatment cost 
##  based on £4400 cost per job (Beatty et al 2011) adjusted to £6193 for Jan 2024
## (from https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator)

## cost per job
cost_job <- 6193

## total cost of intervention
cost_total <- n_treated*cost_job

## cost per qaly gained
cost_qaly <- cost_job/qaly_gain_person

### create summary df
unwtd_df <- data.frame(type="unweighted",
                       measure=c("QALYs gained per person",
                                 "Cost per intervention", 
                                 "ICER"),
                       estimate=c(qaly_gain_person, 
                                  cost_job, cost_qaly))


################################################################################
#####                            PS matched data                           #####
################################################################################
### will need to update - leave out for now #########
#### t0 EQ-5D ------------------------------------------------------------------
#ps_weightit_df_complete <- eq5d_mapper(data=df_matched, 
#                          pcs_var = df_matched$sf12pcs_dv_t0,
#                          mcs_var = df_matched$sf12mcs_dv_t0) %>% 
#  rename(eq5d_t0=eq5d_temp)
#
### create ps data 
#svy_ps_weightit_df_complete <- svydesign(ids = ~1,
#                            data = ps_weightit_df_complete,
#                            weights = ~weights_ps)
#
### create svy table of EQ5D scores by exposure category
#ps_eq5d_t0 <- svyCreateTableOne(vars = "eq5d_t0",
#                                strata = "exposure1",
#                                data = svy_ps_weightit_df_complete,
#                                test = FALSE)

### convert to dataframe with median and IQR
#ps_eq5d_t0 <- data.frame(print(ps_eq5d_t0, nonnormal = "eq5d"))

## check eq-5d total by group
#svyby(~eq5d_t0, ~exposure1, svy_ps_weightit_df_complete, svymean)

### create histogram and boxplots
#par(mfrow=c(2,1))
#svyhist(~eq5d_t0, svy_ps_weightit_df_complete)
#svyboxplot(eq5d_t0~exposure1,svy_ps_weightit_df_complete)

#### t1 EQ-5D ------------------------------------------------------------------
#ps_weightit_df_complete <- eq5d_mapper(data=ps_weightit_df_complete, 
#                          pcs_var = pair_cc_analytic$sf12pcs_dv_t1,
#                          mcs_var = pair_cc_analytic$sf12mcs_dv_t1) %>% 
#  rename(eq5d_t1=eq5d_temp)

### create ps data 
#svy_ps_weightit_df_complete <- svydesign(ids = ~1,
#                            data = ps_weightit_df_complete,
#                            weights = ~weights_ps)

### create svy table of EQ5D scores by exposure category
#ps_eq5d_t1 <- svyCreateTableOne(vars = "eq5d_t1",
#                                strata = "exposure1",
#                                data = svy_ps_weightit_df_complete,
#                                test = FALSE)

### convert to dataframe with median and IQR
#ps_eq5d_t1 <- data.frame(print(ps_eq5d_t1, nonnormal = "eq5d"))

## check eq-5d total by group
#svyby(~eq5d_t1, ~exposure1, svy_ps_weightit_df_complete, svymean)

### create histogram and boxplots
#par(mfrow=c(2,1))
#svyhist(~eq5d_t1, svy_ps_weightit_df_complete)
#svyboxplot(eq5d_t1~exposure1,svy_ps_weightit_df_complete)


#### calculate qalys -----------------------------------------------------------
## this is based on approx 18 months treatment

#ps_qaly_df <- ps_weightit_df_complete %>% mutate(qaly_t0 = 0.5 * eq5d_t0,
#                                    qaly_t1 = 1 * eq5d_t1,
#                                    qaly_total = qaly_t0+qaly_t1)



### create ps data 
#svy_ps_qaly_df <- svydesign(ids = ~1,
#                            data = ps_qaly_df,
#                            weights = ~weights_ps)

### compare qalys for treatment groups
#ps_qaly_grouped <- svyby(~qaly_total, ~exposure1, svy_ps_qaly_df, svytotal)

## add confidence intervals
#ps_qaly_ci <- confint(ps_qaly)
#
#ps_qaly <- ps_qaly %>% bind_cols(ps_qaly_ci) 
#
#ps_qaly <- ps_qaly %>% rename("lci"=`2.5 %`,"uci"=`97.5 %`)

## check qaly totals by group
#svyby(~qaly_t0, ~exposure1, svy_ps_qaly_df, svytotal)
#svyby(~qaly_t1, ~exposure1, svy_ps_qaly_df, svytotal)

## total weighted number of participants
#ps_n_grouped <- svytotal(~exposure1, svy_ps_qaly_df)
#ps_n_grouped <- data.frame(ps_n_grouped) %>% dplyr::select(-SE)

## standardise qalys by sample size of treatment group to get comparable values
#ps_qaly_grouped <- ps_qaly_grouped %>% 
#  bind_cols(ps_n_grouped)


#ps_qaly_grouped <- ps_qaly_grouped %>% mutate(qaly_total_std=qaly_total/total*n_treated)

### calculate QALY gain
#ps_qaly_gain <-  ps_qaly_grouped$qaly_total_std[ps_qaly_grouped$exposure1=="exposed (employed at t1)"] - ps_qaly_grouped$qaly_total_std[ps_qaly_grouped$exposure1=="unexposed"]

### calculate treatment benefit (based on £70k per QALY - UK Govt Green Book)
# https://www.gov.uk/government/publications/the-green-book-appraisal-and-evaluation-in-central-government/the-green-book-2020#valuation-of-costs-and-benefits
#ps_benefit <-  ps_qaly_gain*70000

### calculate treatment cost 
##  based on £4400 cost per job (Beatty et al 2011) adjusted to £6193 for Jan 2024
## (from https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator)

## cost per job
#cost_job <- 6193


## number treated
#ps_n_treated <- ps_n_grouped[2,1]

## total cost of intervention (standardised to unweighted n treated)
#ps_cost_total <- n_treated*cost_job


## cost per qaly gained
#ps_cost_qaly <- ps_cost_total/ps_qaly_gain

### create summary df
#ps_df <- data.frame(type="propensity score matched",
#                    measure=c("QALYs gained", "Treatment benefit", 
#                              "Total cost of intervention", 
#                              "Cost per qaly gained"),
#                    estimate=c(ps_qaly_gain, ps_benefit, 
#                               ps_cost_total, ps_cost_qaly))



################################################################################
#####                                IPTW data                             #####
################################################################################


#### t0 EQ-5D ------------------------------------------------------------------

### create IPTW data 
svy_iptw_weightit_df_complete <- svydesign(ids = ~1,
                              data = weightit_df_complete,
                              weights = ~ps_weights)

### create svy table of EQ5D scores by exposure category
iptw_eq5d_t0 <- svyCreateTableOne(vars = "eq5d_t0",
                                  strata = "exposure1",
                                  data = svy_iptw_weightit_df_complete,
                                  test = FALSE)

### convert to dataframe with median and IQR
iptw_eq5d_t0 <- data.frame(print(iptw_eq5d_t0, nonnormal = "eq5d"))

### create histogram and boxplots
par(mfrow=c(2,1))
svyhist(~eq5d_t0, svy_iptw_weightit_df_complete, breaks = 200)
svyboxplot(eq5d_t0~exposure1,svy_iptw_weightit_df_complete)

#### t1 EQ-5D ------------------------------------------------------------------

### create svy table of EQ5D scores by exposure category
iptw_eq5d_t1 <- svyCreateTableOne(vars = "eq5d_t1",
                                  strata = "exposure1",
                                  data = svy_iptw_weightit_df_complete,
                                  test = FALSE)

### convert to dataframe with median and IQR
iptw_eq5d_t1 <- data.frame(print(iptw_eq5d_t1, nonnormal = "eq5d"))

### create histogram and boxplots
par(mfrow=c(2,1))
svyhist(~eq5d_t1, svy_iptw_weightit_df_complete, breaks = 200)
svyboxplot(eq5d_t1~exposure1,svy_iptw_weightit_df_complete)


#### calculate qalys -----------------------------------------------------------
## this is based on approx 12 months treatment

iptw_qaly_df <- weightit_df_complete %>% mutate(qaly_t0 = eq5d_t0,
                                        qaly_t1 = eq5d_t1,
                                        qaly_diff = (qaly_t1-qaly_t0))



### create iptw data 
svy_iptw_qaly_df <- svydesign(ids = ~1,
                              data = iptw_qaly_df,
                              weights = ~ps_weights)

#### compare qalys for treatment groups
### total difference in QALYs by treatment group
iptw_qaly_grouped <- svyby(~qaly_diff, ~exposure1, svy_iptw_qaly_df, svytotal)
iptw_qaly_grouped <- iptw_qaly_grouped %>% 
  mutate(qaly_diff = qaly_diff/25) %>% dplyr::select(-se)

## add confidence intervals
#iptw_qaly_ci <- confint(iptw_qaly)
#
#iptw_qaly <- iptw_qaly %>% bind_cols(iptw_qaly_ci) 
#
#iptw_qaly <- iptw_qaly %>% rename("lci"=`2.5 %`,"uci"=`97.5 %`)


## total weighted number of participants by treatment group
weighted_n_grouped <- svytotal(~exposure1, svy_iptw_qaly_df)
weighted_n_grouped <- data.frame(weighted_n_grouped) %>% 
  mutate(total = total/25) %>% dplyr::select(-SE)

## calcuate qalys per person
iptw_qaly_grouped <- iptw_qaly_grouped %>% 
  bind_cols(weighted_n_grouped)

## calculate average qalys per person per group to get comparable values
iptw_qaly_grouped <- iptw_qaly_grouped %>% mutate(qaly_person=qaly_diff/total)

### calculate QALY gain
iptw_qaly_gain_person <-  iptw_qaly_grouped$qaly_person[iptw_qaly_grouped$exposure1=="exposed (employed at t1)"] - iptw_qaly_grouped$qaly_person[iptw_qaly_grouped$exposure1=="unexposed"]

### calculate treatment benefit (based on £70k per QALY - UK Govt Green Book)
# https://www.gov.uk/government/publications/the-green-book-appraisal-and-evaluation-in-central-government/the-green-book-2020#valuation-of-costs-and-benefits
iptw_benefit <-  iptw_qaly_gain_person*70000

### calculate treatment cost 
##  based on £4400 cost per job (Beatty et al 2011) adjusted to £6193 for Jan 2024
## (from https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator)

## cost per job
cost_job <- 6193


## number treated
iptw_n_treated <- weighted_n_grouped[2,1]

## total cost of intervention (standardised to unweighted n treated)
iptw_cost_total <- n_treated*cost_job

## calculate average qalys per person per group to get comparable values
iptw_cost_qaly <- cost_job/iptw_qaly_gain_person

### create summary df
iptw_df <- data.frame(type="IPTW",
                      measure=c("QALYs gained per person", 
                                "Cost per intervention", 
                                "ICER"),
                      estimate=c(iptw_qaly_gain_person, 
                                 cost_job, iptw_cost_qaly))


################################################################################
#####                            create qalys df                           #####
################################################################################

qaly_df <- iptw_df %>% 
#  bind_rows(ps_df, df) %>% 
  bind_rows(unwtd_df) %>% 
  pivot_wider(names_from = measure, values_from = estimate) %>% 
  janitor::clean_names() %>% 
  mutate(n_treated = c(iptw_n_treated,
 #                      ps_n_treated,
                       n_treated)) %>% 
  dplyr::select(type, n_treated, everything())

write.csv(qaly_df, "./output/mi/qaly_df.csv")

################################################################################
#####                               scrapbook                              #####
################################################################################

eq5d_min <- −0.148

test <- weightit_df_complete %>% filter(eq5d_t1< eq5d_min)

test2 <- weightit_df_complete %>% filter(pidp=="816041487")
  
