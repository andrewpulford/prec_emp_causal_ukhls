################################################################################

# Precarious employment and health - Understanding Society
# 2-05 - Paired inverse probability of treatment weighted complete case outcome   
# analysis for risk of job loss 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) complete case outcome analysis for IPTW data



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

####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")


#### load IPTW df --------------------------------------------------------------
iptw_df <- readRDS("working_data/cc/weightit_df.rds") # analytic df with IPTW from WeightIt package

iptw_df <- iptw_df %>% 
  dummy_cols(select_columns = c("marital_status_t1",
                                "health_t1"
  ))

iptw_df <- iptw_df %>% janitor::clean_names()


### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
iptw_df$srh_bin_t0 <- factor(iptw_df$srh_bin_t0,
                                levels = c("excellent/very good/good", 
                                           "fair/poor"))
iptw_df$srh_bin_t1 <- factor(iptw_df$srh_bin_t1,
                                levels = c("excellent/very good/good", 
                                           "fair/poor"))


iptw_df$ghq_case4_t0 <- factor(iptw_df$ghq_case4_t0,
                                  levels = c("0-3", "4 or more"))
iptw_df$ghq_case4_t1 <- factor(iptw_df$ghq_case4_t1,
                                  levels = c("0-3", "4 or more"))

iptw_df <- iptw_df %>% 
  mutate(ghq4_outcome_bin = ifelse(ghq_case4_t1=="0-3",0,
                                   ifelse(ghq_case4_t1=="4 or more",1,
                                          "CHECK")))

iptw_df$exposure1 <- factor(iptw_df$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))
iptw_df$exposure2 <- factor(iptw_df$exposure2,
                               levels = c("unexposed",
                                          "exposed (no job loss between t0 and t1"))
# error with exposure 2 to correct at some point

### convert SF-12 outcomes to numeric to allow svyglm to work
iptw_df$sf12pcs_dv_t0 <- as.numeric(iptw_df$sf12pcs_dv_t0)
iptw_df$sf12mcs_dv_t0 <- as.numeric(iptw_df$sf12mcs_dv_t0)
iptw_df$sf12pcs_dv_t1 <- as.numeric(iptw_df$sf12pcs_dv_t1)
iptw_df$sf12mcs_dv_t1 <- as.numeric(iptw_df$sf12mcs_dv_t1)

### ghq-12 change var for checking findings
iptw_df <- iptw_df %>% 
  mutate(ghq_change = paste0(ghq_case4_t0," to ", ghq_case4_t1))

table(iptw_df$ghq_change)

### create weighted data
iptw_svy <- svydesign(ids = ~1,
                         data = iptw_df,
                         weights = ~weightit_ipw)


################################################################################
#####                      weighted regression models                     ######
################################################################################

# initial models prior to doubly robust estimation 

### SF-12 PCS -----------------------
start_time <- Sys.time()
iptw_pcs_glmmTMB_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                      exposure1 +
                                      (1|pidp),
                                    weights = weightit_ipw,
                                    data = iptw_df)

summary(iptw_pcs_glmmTMB_mod)
end_time <- Sys.time()
end_time - start_time

### SF-12 MCS -----------------------
iptw_mcs_glmmTMB_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                      exposure1 +
                                      (1|pidp),
                                 weights = weightit_ipw,
                                 data = iptw_df)

summary(iptw_mcs_glmmTMB_mod)

### poor self-rated health -------------------

iptw_srh_glmmTMB_mod <- glmmTMB(srh_bin_t1 ~ exposure1  +
                                     (1|pidp),
                                   family = binomial(link="logit"),
                                weights = weightit_ipw,
                                data = iptw_df, 
                                   na.action = na.omit)

summary(iptw_srh_glmmTMB_mod)

### GHQ-12 caseness (4+) -------------------

iptw_ghq_glmmTMB_mod <- glmmTMB(ghq_case4_t1 ~ exposure1  +
                                     (1|pidp),
                                   family = binomial(link="logit"),
                                   data = iptw_df, 
                                weights = weightit_ipw,
                                na.action = na.omit)

summary(iptw_ghq_glmmTMB_mod)

## check without multi-level component

ghq_check <- glmmTMB(ghq_case4_t1 ~ exposure1,  
                                family = binomial(link="logit"),
                                data = iptw_df, 
#                                weights = weightit_ipw,
                                na.action = na.omit)

summary(ghq_check)

# add some covariates
ghq_check2 <- glmmTMB(ghq_case4_t1 ~ exposure1 +
                      sex_dv_t0 + age_dv_t0 + degree_bin_t0,  
                     family = binomial(link="logit"),
                     data = iptw_df, 
                     #                                weights = weightit_ipw,
                     na.action = na.omit)

summary(ghq_check2)


################################################################################
#####               double robust weighted regression models               #####
################################################################################



### SF-12 PCS -----------------------
start_time <- Sys.time()
dr_iptw_pcs_glmmTMB_mod <- glmmTMB( sf12pcs_dv_t1 ~
                                      exposure1 +
                                      sex_dv_t0 +
                                         age_dv_t0 +
                                         age_dv_t1 +
                                         non_white_t0 +
                                         marital_status_t0 +
                                         dep_child_bin_t0 +
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
                                         health_t1 +
                                      sf12pcs_dv_t0 +
  # interaction terms
  sex_dv_t0*age_dv_t0 +
  sex_dv_t0*rel_pov_t0 +
  age_dv_t0*rel_pov_t0 +
                                         (1|pidp),
                                    weights = weightit_ipw,
                                    data = iptw_df)
end_time <- Sys.time()
end_time-start_time

dr_iptw_pcs_glmmTMB <- summary(dr_iptw_pcs_glmmTMB_mod)
dr_iptw_pcs_glmmTMB
#fixef(dr_iptw_pcs_glmmTMB_mod)
diagnose(dr_iptw_pcs_glmmTMB_mod)

dr_iptw_pcs_df <- tidy(dr_iptw_pcs_glmmTMB_mod)

## confidence intervals
dr_iptw_pcs_df_ci <- data.frame(confint(dr_iptw_pcs_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_pcs_df_ci <- cbind(rownames(dr_iptw_pcs_df_ci),dr_iptw_pcs_df_ci, row.names=NULL)

dr_iptw_pcs_df_ci <- dr_iptw_pcs_df_ci %>% 
  rename(term = `rownames(dr_iptw_pcs_df_ci)`)

## join dfs together
dr_iptw_pcs_df <- dr_iptw_pcs_df %>% 
  left_join(dr_iptw_pcs_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------
start_time <- Sys.time()
dr_iptw_mcs_glmmTMB_mod <- glmmTMB( sf12mcs_dv_t1 ~
                                      exposure1 +
                                      sex_dv_t0 +
                                      age_dv_t0 +
                                      age_dv_t1 +
                                      non_white_t0 +
                                      marital_status_t0 +
                                      dep_child_bin_t0 +
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
                                      health_t1 +
                                      sf12mcs_dv_t0 +
  # interaction terms
  sex_dv_t0*age_dv_t0 +
  sex_dv_t0*rel_pov_t0 +
  age_dv_t0*rel_pov_t0 +
                                         (1|pidp),
                                    weights = weightit_ipw,
                                    data = iptw_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_mcs_glmmTMB_mod)

diagnose(dr_iptw_mcs_glmmTMB_mod)

dr_iptw_mcs_df <- tidy(dr_iptw_mcs_glmmTMB_mod)

## confidence intervals
dr_iptw_mcs_df_ci <- data.frame(confint(dr_iptw_mcs_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_mcs_df_ci <- cbind(rownames(dr_iptw_mcs_df_ci),dr_iptw_mcs_df_ci, row.names=NULL)

dr_iptw_mcs_df_ci <- dr_iptw_mcs_df_ci %>% 
  rename(term = `rownames(dr_iptw_mcs_df_ci)`)

## join dfs together
dr_iptw_mcs_df <- dr_iptw_mcs_df %>% 
  left_join(dr_iptw_mcs_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

start_time <- Sys.time()
dr_iptw_srh_glmmTMB_mod <- glmmTMB( srh_bin_t1 ~
                                      exposure1 +
                                      sex_dv_t0 +
                                      age_dv_t0 +
                                      age_dv_t1 +
                                      non_white_t0 +
                                      marital_status_t0 +
                                      dep_child_bin_t0 +
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
                                      health_t1 +
                                      srh_bin_t0 +
  # interaction terms
  sex_dv_t0*age_dv_t0 +
  sex_dv_t0*rel_pov_t0 +
  age_dv_t0*rel_pov_t0 +
                                         (1|pidp),
                                       family = binomial(link="logit"),
                                    weights = weightit_ipw,
                                    data = iptw_df)
end_time <- Sys.time()
end_time-start_time

summary(dr_iptw_srh_glmmTMB_mod)
diagnose(dr_iptw_srh_glmmTMB_mod)

dr_iptw_srh_df <- tidy(dr_iptw_srh_glmmTMB_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_srh_df_ci <- data.frame(confint(dr_iptw_srh_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_srh_df_ci <- cbind(rownames(dr_iptw_srh_df_ci),dr_iptw_srh_df_ci, row.names=NULL)

dr_iptw_srh_df_ci <- dr_iptw_srh_df_ci %>% 
  rename(term = `rownames(dr_iptw_srh_df_ci)`)

## join dfs together
dr_iptw_srh_df <- dr_iptw_srh_df %>% 
  left_join(dr_iptw_srh_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------
start_time <- Sys.time()
dr_iptw_ghq_glmmTMB_mod <- glmmTMB( ghq_case4_t1 ~
                                      exposure1 +
                                      sex_dv_t0 +
                                      age_dv_t0 +
                                      age_dv_t1 +
                                      non_white_t0 +
                                      marital_status_t0 +
                                      dep_child_bin_t0 +
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
                                      health_t1 +
                                      ghq_case4_t0 +
  # interaction terms
  sex_dv_t0*age_dv_t0 +
  sex_dv_t0*rel_pov_t0 +
  age_dv_t0*rel_pov_t0 +
                                           (1|pidp),
                                       family = binomial(link="logit"),
                                    weights = weightit_ipw,
                                    data = iptw_df)
end_time <- Sys.time()
end_time-start_time

dr_iptw_ghq_glmmTMB_summary <- summary(dr_iptw_ghq_glmmTMB_mod)

write_rds(dr_iptw_ghq_glmmTMB_summary, "./output/scrapbook/dr_iptw_ghq_glmmTMB_summary.rds")

diagnose(dr_iptw_ghq_glmmTMB_mod)

dr_iptw_ghq_df <- tidy(dr_iptw_ghq_glmmTMB_mod) %>% 
  mutate(estimate = exp(estimate))


## confidence intervals
dr_iptw_ghq_df_ci <- data.frame(confint(dr_iptw_ghq_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)  %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_iptw_ghq_df_ci <- cbind(rownames(dr_iptw_ghq_df_ci),dr_iptw_ghq_df_ci, row.names=NULL)

dr_iptw_ghq_df_ci <- dr_iptw_ghq_df_ci %>% 
  rename(term = `rownames(dr_iptw_ghq_df_ci)`)

## join dfs together
dr_iptw_ghq_df <- dr_iptw_ghq_df %>% 
  left_join(dr_iptw_ghq_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### create single summary df for iptw double-robust cc outcomes ------------

dr_iptw_df <- dr_iptw_pcs_df %>% 
  bind_rows(dr_iptw_mcs_df,
            dr_iptw_srh_df,
            dr_iptw_ghq_df) %>% 
  filter(term=="exposed (employed at t1)") %>% 
  dplyr::select(-c(term, Estimate, group, component)) %>% 
  dplyr::select(outcome, effect, est_type, estimate, std.error, p.value, lci, uci)

write.csv(dr_iptw_df, "./output/cc/weighted_outcomes/dr_iptw_df.csv")


diagnose(dr_iptw_pcs_glmmTMB_mod)
diagnose(dr_iptw_mcs_glmmTMB_mod)
diagnose(dr_iptw_srh_glmmTMB_mod)
diagnose(dr_iptw_ghq_glmmTMB_mod)

############## GHQ-12 checks


## descriptive table
iptw_df %>% group_by(exposure1,ghq_case4_t0,ghq_case4_t1) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(exposure1) %>% 
  mutate(d=sum(n)) %>% 
  ungroup() %>% 
  mutate(pc=n/d*100)

## descriptive table2
iptw_df %>% group_by(exposure1,ghq_case4_t1) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(exposure1) %>% 
  mutate(d=sum(n)) %>% 
  ungroup() %>% 
  mutate(pc=n/d*100)

## manual OR calculation
odds_exp <- 14724/74356
odd_unexp <- 605/1001
OR <- odds_exp/odd_unexp

## table one
ghq_change_summary <- svyCreateTableOne(vars=c("ghq_case4_t0","ghq_case4_t1","ghq_change"),
                  strata="exposure1",
                  data=iptw_svy)

ghq_change_summary_sav <- print(ghq_change_summary)

write.csv(ghq_change_summary_sav, "./output/scrapbook/ghq_change_summary_cc.csv")

##### regression model for GHQ-12 scale
start_time <- Sys.time()
dr_iptw_ghq_scale_glmmTMB_mod <- glmmTMB( scghq2_dv_t1 ~
                                            exposure1 +
                                            sex_dv_t0 +
                                            age_dv_t0 +
                                            age_dv_t1 +
                                            non_white_t0 +
                                            marital_status_t0 +
                                            dep_child_bin_t0 +
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
                                            health_t1 +
                                            scghq2_dv_t0 +
                                      # interaction terms
                                      sex_dv_t0*age_dv_t0 +
                                      sex_dv_t0*rel_pov_t0 +
                                      age_dv_t0*rel_pov_t0 +
                                      (1|pidp),
                                    weights = weightit_ipw,
                                    data = iptw_df)
end_time <- Sys.time()
end_time-start_time

dr_iptw_ghq_scale_glmmTMB <- summary(dr_iptw_ghq_scale_glmmTMB_mod)
dr_iptw_ghq_scale_glmmTMB
#fixef(dr_iptw_ghq_scale_glmmTMB_mod)
diagnose(dr_iptw_ghq_scale_glmmTMB_mod)

dr_iptw_ghq_scale_df <- tidy(dr_iptw_ghq_scale_glmmTMB_mod)

## confidence intervals
dr_iptw_ghq_scale_df_ci <- data.frame(confint(dr_iptw_ghq_scale_glmmTMB_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
dr_iptw_ghq_scale_df_ci <- cbind(rownames(dr_iptw_ghq_scale_df_ci),dr_iptw_ghq_scale_df_ci, row.names=NULL)

dr_iptw_ghq_scale_df_ci <- dr_iptw_ghq_scale_df_ci %>% 
  rename(term = `rownames(dr_iptw_ghq_scale_df_ci)`)

## join dfs together
dr_iptw_ghq_scale_df <- dr_iptw_ghq_scale_df %>% 
  left_join(dr_iptw_ghq_scale_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 scale",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

##### basic logistic regression models (no multi-level component)

glm(ghq_case4_t1 ~ exposure1 , family = binomial(link = "logit"), 
    data = iptw_df)
exp(-1.1159)

################################################################################
#####                    GHQ-12 sub-group analyses                        ######
################################################################################

#### create df's ---------------------------------------------------------------
female <- iptw_df %>% filter(sex_dv_t0_female==1)
male <- iptw_df %>% filter(sex_dv_t0_male==1)

#### sex -----------------------------------------------------------------------
### female
start_time <- Sys.time()
dr_iptw_ghq_glmmTMB_mod <- glmmTMB( ghq_case4_t1 ~
                                      exposure1 +
 #                                     sex_dv_t0 +
                                      age_dv_t0 +
                                      age_dv_t1 +
                                      non_white_t0 +
                                      marital_status_t0 +
                                      dep_child_bin_t0 +
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
                                      health_t1 +
                                      ghq_case4_t0 +
                                      # interaction terms
 #                                     sex_dv_t0*age_dv_t0 +
 #                                     sex_dv_t0*rel_pov_t0 +
                                      age_dv_t0*rel_pov_t0 +
                                      (1|pidp),
                                    family = binomial(link="logit"),
                                    weights = weightit_ipw,
                                    data = female)
end_time <- Sys.time()
end_time-start_time

dr_iptw_ghq_glmmTMB_summary <- summary(dr_iptw_ghq_glmmTMB_mod)


