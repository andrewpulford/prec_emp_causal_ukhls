################################################################################

# Precarious employment and health - Understanding Society
# 2-06 - complete case EQ5D/QALY analysis   
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
library(fastDummies) # for creating dummy variables
library(Matching) # PS matching
library(tableone) # for creating table one
library(survey) # for PS weighting
library(reshape2) # Reorganizing data
library(broom) # for tidying regression outputs into df format
library(broom.mixed) # for tidying MLM regression outputs into df format
library(glmmTMB) # for multi-level modelling (faster than lme4)
#library(doParallel)
#library(foreach)
#cluster = makeCluster(detectCores() - 1) # convention to leave 1 core for OS
#registerDoParallel(cluster)

################################################################################
#####                         load and prepare data                        #####
################################################################################


#### load unweighted analytic df -----------------------------------------------
pair_cc_analytic <- readRDS("./working_data/cc/pair_cc_analytic.rds")

### convert SF-12 outcomes to numeric
pair_cc_analytic$sf12pcs_dv_t0 <- as.character(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t0 <- as.character(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t1 <- as.character(pair_cc_analytic$sf12pcs_dv_t1)
pair_cc_analytic$sf12mcs_dv_t1 <- as.character(pair_cc_analytic$sf12mcs_dv_t1)

pair_cc_analytic$sf12pcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t1)
pair_cc_analytic$sf12mcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t1)


#### load PR matched analytic df -----------------------------------------------
df_matched <- readRDS("working_data/cc/matchit_df.rds") # PS matched data

### convert SF-12 outcomes to numeric
df_matched$sf12pcs_dv_t0 <- as.character(df_matched$sf12pcs_dv_t0)
df_matched$sf12mcs_dv_t0 <- as.character(df_matched$sf12mcs_dv_t0)
df_matched$sf12pcs_dv_t1 <- as.character(df_matched$sf12pcs_dv_t1)
df_matched$sf12mcs_dv_t1 <- as.character(df_matched$sf12mcs_dv_t1)

df_matched$sf12pcs_dv_t0 <- as.numeric(df_matched$sf12pcs_dv_t0)
df_matched$sf12mcs_dv_t0 <- as.numeric(df_matched$sf12mcs_dv_t0)
df_matched$sf12pcs_dv_t1 <- as.numeric(df_matched$sf12pcs_dv_t1)
df_matched$sf12mcs_dv_t1 <- as.numeric(df_matched$sf12mcs_dv_t1)

#### load IPTW analytic df -----------------------------------------------------
iptw_df <- readRDS("working_data/cc/weightit_df.rds") # analytic df with IPTW from WeightIt package

### convert SF-12 outcomes to numeric
iptw_df$sf12pcs_dv_t0 <- as.character(iptw_df$sf12pcs_dv_t0)
iptw_df$sf12mcs_dv_t0 <- as.character(iptw_df$sf12mcs_dv_t0)
iptw_df$sf12pcs_dv_t1 <- as.character(iptw_df$sf12pcs_dv_t1)
iptw_df$sf12mcs_dv_t1 <- as.character(iptw_df$sf12mcs_dv_t1)

iptw_df$sf12pcs_dv_t0 <- as.numeric(iptw_df$sf12pcs_dv_t0)
iptw_df$sf12mcs_dv_t0 <- as.numeric(iptw_df$sf12mcs_dv_t0)
iptw_df$sf12pcs_dv_t1 <- as.numeric(iptw_df$sf12pcs_dv_t1)
iptw_df$sf12mcs_dv_t1 <- as.numeric(iptw_df$sf12mcs_dv_t1)


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
### call function
unwtd_eq5d_df <- eq5d_mapper(data=pair_cc_analytic, 
                             pcs_var = pair_cc_analytic$sf12pcs_dv_t0,
                             mcs_var = pair_cc_analytic$sf12mcs_dv_t0) %>% 
  rename(eq5d_t0=eq5d_temp)

### summary stats
summary(unwtd_eq5d_df$eq5d_t0)

### check distribution with histogram
hist(unwtd_eq5d_df$eq5d_t0, breaks = 200)

### plot E8-5D scores against S-12 PCS and MCS
plot(unwtd_eq5d_df$sf12pcs_dv_t0,unwtd_eq5d_df$eq5d_t0)
plot(unwtd_eq5d_df$sf12mcs_dv_t0,unwtd_eq5d_df$eq5d_t0)

### calculate mean EQ-5D by treatment group
unwtd_eq5d_df %>% group_by(exposure1) %>% summarise(mean(eq5d_t0))

### EQ-5D distribution by treatment group
boxplot(unwtd_eq5d_df$eq5d_t0 ~ unwtd_eq5d_df$exposure1)

#### t1 EQ-5D ------------------------------------------------------------------
### call function
unwtd_eq5d_df <- eq5d_mapper(data=unwtd_eq5d_df, 
                             pcs_var = pair_cc_analytic$sf12pcs_dv_t1,
                             mcs_var = pair_cc_analytic$sf12mcs_dv_t1) %>% 
  rename(eq5d_t1=eq5d_temp)

### summary stats
summary(unwtd_eq5d_df$eq5d_t1)

### check distribution with histogram
hist(unwtd_eq5d_df$eq5d_t1, breaks = 200)

### plot E8-5D scores against S-12 PCS and MCS
plot(unwtd_eq5d_df$sf12pcs_dv_t1,unwtd_eq5d_df$eq5d_t1)
plot(unwtd_eq5d_df$sf12mcs_dv_t1,unwtd_eq5d_df$eq5d_t1)

### calculate mean EQ-5D by treatment group
unwtd_eq5d_df %>% group_by(exposure1) %>% summarise(mean(eq5d_t1))

### EQ-5D distribution by treatment group
boxplot(unwtd_eq5d_df$eq5d_t1 ~ unwtd_eq5d_df$exposure1)


#### calculate qalys -----------------------------------------------------------
## this is based on approx 12 months treatment
## QALY at single time-point is equal to EQ-5D score at same time-point
## for analysis we calculate the difference between QALY at t1 and t0
unwtd_qaly_df <- unwtd_eq5d_df %>% 
  mutate(qaly_t0 = eq5d_t0,
         qaly_t1 = eq5d_t1,
         qaly_diff = qaly_t1-qaly_t0) 

#### compare qalys for treatment groups
### total difference in QALYs by treatment group
unwtd_qaly_grouped <- unwtd_qaly_df %>% group_by(exposure1) %>% 
  summarise(qaly_diff_total = sum(qaly_diff))

### total number of participants with treatment group
unwtd_n_grouped <- unwtd_qaly_df %>% group_by(exposure1) %>% 
  summarise(n = n())

## total number in treatment group
unwtd_n_treated <- sum(unwtd_qaly_df$exposure1=="exposed (employed at t1)")

## calculate average qalys per person per group to get comparable values
unwtd_qaly_grouped <- unwtd_qaly_grouped %>% 
  left_join(unwtd_n_grouped) %>% 
  mutate(qaly_diff_person=qaly_diff_total/n)


### calculate QALY gain
## QALY gain per treated person
unwtd_qaly_gain_person <- unwtd_qaly_grouped$qaly_diff_person[unwtd_qaly_grouped$exposure1=="exposed (employed at t1)"] - unwtd_qaly_grouped$qaly_diff_person[unwtd_qaly_grouped$exposure1=="unexposed"]

## total QALYs gained for treatment group
qualy_gain_total <- unwtd_qaly_gain_person*unwtd_n_treated

### calculate treatment benefit (based on £70k per QALY - UK Govt Green Book)
# https://www.gov.uk/government/publications/the-green-book-appraisal-and-evaluation-in-central-government/the-green-book-2020#valuation-of-costs-and-benefits
unwtd_benefit <-  unwtd_qaly_gain_person*70000

### calculate treatment cost 
##  based on £4400 cost per job (Beatty et al 2011) adjusted to £6193 for Jan 2024
## (from https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator)

## cost per job
cost_job <- 6193

## total cost of intervention
unwtd_cost_total <- unwtd_n_treated*cost_job

## cost per qaly gained (ICER)
unwtd_icer <- cost_job/unwtd_qaly_gain_person

### create summary df
unwtd_df <- data.frame(type="unweighted",
                       measure=c("QALYs gained per person", 
                                 "Cost per intervention", 
                                 "ICER"),
                       estimate=c(unwtd_qaly_gain_person, 
                                  cost_job, unwtd_icer))


################################################################################
#####                            PS matched data                           #####
################################################################################

#### t0 EQ-5D ------------------------------------------------------------------
ps_eq5d_df <- eq5d_mapper(data=df_matched, 
                          pcs_var = df_matched$sf12pcs_dv_t0,
                          mcs_var = df_matched$sf12mcs_dv_t0) %>% 
  rename(eq5d_t0=eq5d_temp)

### create ps data 
svy_ps_eq5d_df <- svydesign(ids = ~1,
                            data = ps_eq5d_df,
                            weights = ~weights_ps)

### create svy table of EQ5D scores by exposure category
ps_eq5d_t0 <- svyCreateTableOne(vars = "eq5d_t0",
                                strata = "exposure1",
                                data = svy_ps_eq5d_df,
                                test = FALSE)

### convert to dataframe with median and IQR
ps_eq5d_t0 <- data.frame(print(ps_eq5d_t0, nonnormal = "eq5d"))

## check eq-5d total by group
svyby(~eq5d_t0, ~exposure1, svy_ps_eq5d_df, svymean)

### create histogram and boxplots
par(mfrow=c(2,1))
svyhist(~eq5d_t0, svy_ps_eq5d_df, breaks = 200)
svyboxplot(eq5d_t0~exposure1,svy_ps_eq5d_df)

#### t1 EQ-5D ------------------------------------------------------------------
ps_eq5d_df <- eq5d_mapper(data=ps_eq5d_df, 
                          pcs_var = pair_cc_analytic$sf12pcs_dv_t1,
                          mcs_var = pair_cc_analytic$sf12mcs_dv_t1) %>% 
  rename(eq5d_t1=eq5d_temp)

### create ps data 
svy_ps_eq5d_df <- svydesign(ids = ~1,
                            data = ps_eq5d_df,
                            weights = ~weights_ps)

### create svy table of EQ5D scores by exposure category
ps_eq5d_t1 <- svyCreateTableOne(vars = "eq5d_t1",
                                strata = "exposure1",
                                data = svy_ps_eq5d_df,
                                test = FALSE)

### convert to dataframe with median and IQR
ps_eq5d_t1 <- data.frame(print(ps_eq5d_t1, nonnormal = "eq5d"))

## check eq-5d total by group
svyby(~eq5d_t1, ~exposure1, svy_ps_eq5d_df, svymean)

### create histogram and boxplots
par(mfrow=c(2,1))
svyhist(~eq5d_t1, svy_ps_eq5d_df, breaks = 200)
svyboxplot(eq5d_t1~exposure1,svy_ps_eq5d_df)


#### calculate qalys -----------------------------------------------------------
## this is based on approx 12 months treatment

ps_qaly_df <- ps_eq5d_df %>% mutate(qaly_t0 = eq5d_t0,
                                    qaly_t1 = eq5d_t1,
                                    qaly_diff = qaly_t1-qaly_t0)



### create ps data 
svy_ps_qaly_df <- svydesign(ids = ~1,
                            data = ps_qaly_df,
                            weights = ~weights_ps)


#### compare qalys for treatment groups
### total difference in QALYs by treatment group
ps_qaly_grouped <- svyby(~qaly_diff, ~exposure1, svy_ps_qaly_df, svytotal)


## total weighted number of participants by treatment group
ps_n_grouped <- svytotal(~exposure1, svy_ps_qaly_df)
ps_n_grouped <- data.frame(ps_n_grouped) %>% dplyr::select(-SE)

## calcuate qalys per person
ps_qaly_grouped <- ps_qaly_grouped %>% 
  bind_cols(ps_n_grouped)

## calculate average qalys per person per group to get comparable values
ps_qaly_grouped <- ps_qaly_grouped %>% mutate(qaly_person=qaly_diff/total)

### calculate QALY gain
## QALY gain per treated person
ps_qaly_gain_person <-  ps_qaly_grouped$qaly_person[ps_qaly_grouped$exposure1=="exposed (employed at t1)"] - ps_qaly_grouped$qaly_person[ps_qaly_grouped$exposure1=="unexposed"]

### calculate treatment benefit (based on £70k per QALY - UK Govt Green Book)
# https://www.gov.uk/government/publications/the-green-book-appraisal-and-evaluation-in-central-government/the-green-book-2020#valuation-of-costs-and-benefits
ps_benefit <-  ps_qaly_gain_person*70000

### calculate treatment cost 
##  based on £4400 cost per job (Beatty et al 2011) adjusted to £6193 for Jan 2024
## (from https://www.bankofengland.co.uk/monetary-policy/inflation/inflation-calculator)

## cost per job
cost_job <- 6193


## number treated
ps_n_treated <- ps_n_grouped[2,1]

## total cost of intervention 
ps_cost_total <- ps_n_treated*cost_job


## cost per qaly gained (ICER)
ps_cost_qaly <- cost_job/ps_qaly_gain_person

### create summary df
ps_df <- data.frame(type="propensity score matched",
                    measure=c("QALYs gained per person", 
                              "Cost per intervention", 
                              "ICER"),
                    estimate=c(ps_qaly_gain_person, 
                               cost_job, ps_cost_qaly))



################################################################################
#####                                IPTW data                             #####
################################################################################

#### t0 EQ-5D ------------------------------------------------------------------
iptw_eq5d_df <- eq5d_mapper(data=iptw_df, 
                            pcs_var = iptw_df$sf12pcs_dv_t0,
                            mcs_var = iptw_df$sf12mcs_dv_t0) %>% 
  rename(eq5d_t0=eq5d_temp)

### create IPTW data 
svy_iptw_eq5d_df <- svydesign(ids = ~1,
                              data = iptw_eq5d_df,
                              weights = ~weightit_ipw)

### create svy table of EQ5D scores by exposure category
iptw_eq5d_t0 <- svyCreateTableOne(vars = "eq5d_t0",
                                  strata = "exposure1",
                                  data = svy_iptw_eq5d_df,
                                  test = FALSE)

### convert to dataframe with median and IQR
iptw_eq5d_t0 <- data.frame(print(iptw_eq5d_t0, nonnormal = "eq5d"))

### create histogram and boxplots
par(mfrow=c(2,1))
svyhist(~eq5d_t0, svy_iptw_eq5d_df, breaks = 200)
svyboxplot(eq5d_t0~exposure1,svy_iptw_eq5d_df)

#### t1 EQ-5D ------------------------------------------------------------------
iptw_eq5d_df <- eq5d_mapper(data=iptw_eq5d_df, 
                            pcs_var = pair_cc_analytic$sf12pcs_dv_t1,
                            mcs_var = pair_cc_analytic$sf12mcs_dv_t1) %>% 
  rename(eq5d_t1=eq5d_temp)

### create IPTW data 
svy_iptw_eq5d_df <- svydesign(ids = ~1,
                              data = iptw_eq5d_df,
                              weights = ~weightit_ipw)

### create svy table of EQ5D scores by exposure category
iptw_eq5d_t1 <- svyCreateTableOne(vars = "eq5d_t1",
                                  strata = "exposure1",
                                  data = svy_iptw_eq5d_df,
                                  test = FALSE)

### convert to dataframe with median and IQR
iptw_eq5d_t1 <- data.frame(print(iptw_eq5d_t1, nonnormal = "eq5d"))

### create histogram and boxplots
par(mfrow=c(2,1))
svyhist(~eq5d_t1, svy_iptw_eq5d_df, breaks = 200)
svyboxplot(eq5d_t1~exposure1,svy_iptw_eq5d_df)


#### calculate qalys -----------------------------------------------------------
## this is based on approx 12 months treatment

iptw_qaly_df <- iptw_eq5d_df %>% mutate(qaly_t0 = eq5d_t0,
                                        qaly_t1 = eq5d_t1,
                                        qaly_diff = qaly_t1-qaly_t0)



### create iptw data 
svy_iptw_qaly_df <- svydesign(ids = ~1,
                              data = iptw_qaly_df,
                              weights = ~weightit_ipw)

#### compare qalys for treatment groups
### total difference in QALYs by treatment group
iptw_qaly_grouped <- svyby(~qaly_diff, ~exposure1, svy_iptw_qaly_df, svytotal)

## add confidence intervals
#iptw_qaly_ci <- confint(iptw_qaly)
#
#iptw_qaly <- iptw_qaly %>% bind_cols(iptw_qaly_ci) 
#
#iptw_qaly <- iptw_qaly %>% rename("lci"=`2.5 %`,"uci"=`97.5 %`)

## total weighted number of participants by treatment group
weighted_n_grouped <- svytotal(~exposure1, svy_iptw_qaly_df)
weighted_n_grouped <- data.frame(weighted_n_grouped) %>% dplyr::select(-SE)

## calcuate qalys per person
iptw_qaly_grouped <- iptw_qaly_grouped %>% 
  bind_cols(weighted_n_grouped)

## calculate average qalys per person per group to get comparable values
iptw_qaly_grouped <- iptw_qaly_grouped %>% mutate(qaly_person=qaly_diff/total)

### calculate QALY gain
## QALY gain per treated person
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

## total cost of intervention
iptw_cost_total <- iptw_n_treated*cost_job

## cost per qaly gained (ICER)
iptw_cost_qaly <- cost_job/iptw_qaly_gain_person

### create summary df
iptw_df <- data.frame(type="IPTW",
                      measure=c("QALYs gained per person", 
                                "Cost per intervention", 
                                "ICER"),
                      estimate=c(iptw_qaly_gain_person, 
                                 cost_job, iptw_cost_qaly))



################################################################################
##### create qalys df #####
################################################################################

qaly_df <- iptw_df %>% 
  bind_rows(ps_df, unwtd_df) %>% 
  pivot_wider(names_from = measure, values_from = estimate) %>% 
  janitor::clean_names() %>% 
  mutate(n_treated = c(iptw_n_treated,ps_n_treated,unwtd_n_treated)) %>% 
  dplyr::select(type, n_treated, everything())

write.csv(qaly_df, "./output/cc/qaly_df.csv")

################################################################################
##### scrapbook ######
################################################################################

#### unweighted ----------------------------------------------------------------

## Franks
eq5d_test <- pair_cc_analytic %>% 
  mutate(eq5d_t0 = 0.84690 + 
         (sf12pcs_dv_t0-49.9)*0.01261 + 
         (sf12mcs_dv_t0-51.5)*0.00759 -
         (sf12pcs_dv_t0-49.9)^2*0.00009 - 
         (sf12mcs_dv_t0-51.5)^2*0.00015 - 
         (sf12pcs_dv_t0-49.9)*(sf12mcs_dv_t0-51.5)*0.00015) %>% 
  mutate(eq5d_t1 = 0.84690 + 
           (sf12pcs_dv_t1-49.9)*0.01261 + 
           (sf12mcs_dv_t1-51.5)*0.00759 -
           (sf12pcs_dv_t1-49.9)^2*0.00009 - 
           (sf12mcs_dv_t1-51.5)^2*0.00015 - 
           (sf12pcs_dv_t1-49.9)*(sf12mcs_dv_t1-51.5)*0.00015)

## Lawrence
lawrence_test <- pair_cc_analytic %>% 
  mutate(eq5d_t0 = -1.6984 +
           (sf12pcs_dv_t0)*0.07927 +
           (sf12mcs_dv_t0)*0.02859 +
           (sf12pcs_dv_t0)*(sf12mcs_dv_t0)*-0.000126 +
           (sf12pcs_dv_t0)^2*-0.00141 +
           (sf12mcs_dv_t0)^2*-0.00014 +
           (sf12pcs_dv_t0)^3*0.0000107) %>% 
  mutate(eq5d_t1 = -1.6984 +
           (sf12pcs_dv_t1)*0.07927 +
           (sf12mcs_dv_t1)*0.02859 +
           (sf12pcs_dv_t1)*(sf12mcs_dv_t0)*-0.000126 +
           (sf12pcs_dv_t1)^2*-0.00141 +
           (sf12mcs_dv_t1)^2*-0.00014 +
           (sf12pcs_dv_t1)^3*0.0000107) 


summary(lawrence_test$eq5d_t0)
hist(lawrence_test$eq5d_t0, breaks = 200)


plot(lawrence_test$eq5d_t0,lawrence_test$sf12pcs_dv_t0)

lawrence_test %>% group_by(exposure1) %>% summarise(mean(eq5d_t0))

boxplot(lawrence_test$eq5d_t0 ~ lawrence_test$exposure1)

# simpler Lawrence
lawrence_test2 <- pair_cc_analytic %>% 
  mutate(eq5d_t0 = -0.372 +
           (sf12pcs_dv_t0)*0.01411 +
           (sf12mcs_dv_t0)*0.00967) %>%  
  mutate(eq5d_t1 = -0.372 +
           (sf12pcs_dv_t1)*0.01411 +
           (sf12mcs_dv_t1)*0.00967)

summary(lawrence_test2$eq5d_t0)
hist(lawrence_test2$eq5d_t0, breaks = 200)


plot(lawrence_test2$eq5d_t0,lawrence_test2$sf12pcs_dv_t0)

lawrence_test2 %>% group_by(exposure1) %>% summarise(mean(eq5d_t0))

boxplot(lawrence_test2$eq5d_t0 ~ lawrence_test2$exposure1)

### qaly calculations
## Franks 
#1.	Calculate total QALYs at T0 by exposure group
#2.	Calculate total QALYs at T1 by exposure group
unwtd_qaly_df2 <- eq5d_test %>% 
  group_by(exposure1) %>% 
  summarise(total_qaly_t0 = sum(eq5d_t0),
            total_qaly_t1 = sum(eq5d_t1),
            d = n()) %>% 
#3.	Deduct 1 from 2 to give total QALY change by exposure group
  mutate(total_qaly_diff = total_qaly_t1-total_qaly_t0,
#4.	Divide 3 by the group size to give QALY change per person by exposure group
         qaly_diff_person = total_qaly_diff/d)
#5.	Compare 4 between exposure groups to give average treatment effect. 
ate <- unwtd_qaly_df2$qaly_diff_person[unwtd_qaly_df2$exposure1=="exposed (employed at t1)"] - unwtd_qaly_df2$qaly_diff_person[unwtd_qaly_df2$exposure1=="unexposed"]

## Lawrence 
#1.	Calculate total QALYs at T0 by exposure group
#2.	Calculate total QALYs at T1 by exposure group
unwtd_qaly_df3 <- lawrence_test %>% 
  group_by(exposure1) %>% 
  summarise(total_qaly_t0 = sum(eq5d_t0),
            total_qaly_t1 = sum(eq5d_t1),
            d = n()) %>% 
  #3.	Deduct 1 from 2 to give total QALY change by exposure group
  mutate(total_qaly_diff = total_qaly_t1-total_qaly_t0,
         #4.	Divide 3 by the group size to give QALY change per person by exposure group
         qaly_diff_person = total_qaly_diff/d)
#5.	Compare 4 between exposure groups to give average treatment effect. 
ate2 <- unwtd_qaly_df3$qaly_diff_person[unwtd_qaly_df3$exposure1=="exposed (employed at t1)"] - unwtd_qaly_df3$qaly_diff_person[unwtd_qaly_df3$exposure1=="unexposed"]

max(lawrence_test$eq5d_t0)
max(lawrence_test$eq5d_t1)

lawrence_test2 %>% filter(eq5d_t0>1) %>% 
  dplyr::select(sf12pcs_dv_t0,sf12mcs_dv_t0,eq5d_t0)


#### checking against MI dataset -----------------------------------------------

test2 <- pair_cc_analytic %>% filter(pidp=="816041487")

