################################################################################

# Precarious employment and health - Understanding Society
# 3-05 - Paired propensity weighted complete case outcome analysis for risk of 
# job loss 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) 



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
library(Matching) # PS matching
library(tableone) # for creating table one
library(survey) # for PS weighting
library(reshape2) # Reorganizing data
library(broom) # for tidying regression outputs into df format
#library(doParallel)
#library(foreach)
#cluster = makeCluster(detectCores() - 1) # convention to leave 1 core for OS
#registerDoParallel(cluster)

################################################################################
#####                         load and prepare data                        #####
################################################################################

#### create vectors for variables ----------------------------------------------
## id and weights vars
id_wt_vector <- c("pidp", "psu", "strata", "wt_name", "wt_value")

## baseline covariates
cov_vector <- c("sex_dv_t0", 
                "age_dv_t0",  
                "non_white_t0", 
                "marital_status_t0",
                "hiqual_dv_t0", 
                "gor_dv_t0",
                "sic2007_section_lab_t0",
                "soc2000_major_group_title_t0",
                "jbft_dv_t0",
                "small_firm_t0",
                "emp_contract_t0",
                "broken_emp_t0",
                "j2has_dv_t0",
                "rel_pov_t0",
#                "jbhrs_t0",
                "health_t0",
                "srh_bin_t0",
                "ghq_case4_t0",
                "sf12mcs_dv_t0",
                "sf12pcs_dv_t0")
# use missing cat
# don't create separate outcome df's unless high # missing

## time varying covariates at t1
cov_vector2 <- c("age_dv_t1",  
                 "marital_status_t1",
                 "gor_dv_t1",
                 "health_t1",
                 "srh_bin_t1",
                 "ghq_case4_t1",
                 "sf12mcs_dv_t1",
                 "sf12pcs_dv_t1")

## outcome vector including exposure vars
outcome_vector <- c("srh_bin_t1",
                    "ghq_case4_t1",
                    "sf12mcs_dv_t1",
                    "sf12pcs_dv_t1",
                    "exposure1",
                    "exposure2")


outcome_vector2 <- c("sf12pcs_dv_t1",
                     "sf12mcs_dv_t1",
                     "srh_bin_t1",
                     "ghq_case4_t1")

#### load propsenty matching df for weights ------------------------------------
pair_cc_ps <- readRDS("working_data/pair_cc_ps.rds")

### create weighted data
pair_cc_ps_svy <- svydesign(ids = ~1,
                            data = pair_cc_ps,
                            weights = ~mw)



####load eligible cases --------------------------------------------------------
#pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds") %>% 
#  dplyr::select(c(id_wt_vector, cov_vector, cov_vector2, outcome_vector)) %>% 
#  dplyr::select(-c(psu, strata, wt_name, wt_value))

### convert binary outcome and exposure vars to factors and relevel to allow svyglm to work
pair_cc_ps$srh_bin_t0 <- factor(pair_cc_ps$srh_bin_t0,
                                      levels = c("good/fair/poor", 
                                                 "excellent/very good"))
pair_cc_ps$srh_bin_t1 <- factor(pair_cc_ps$srh_bin_t1,
                                      levels = c("good/fair/poor", 
                                                 "excellent/very good"))

pair_cc_ps$ghq_case4_t0 <- factor(pair_cc_ps$ghq_case4_t0,
                                        levels = c("0-3", "4 or more"))
pair_cc_ps$ghq_case4_t1 <- factor(pair_cc_ps$ghq_case4_t1,
                                        levels = c("0-3", "4 or more"))


pair_cc_ps$exposure1 <- factor(pair_cc_ps$exposure1,
                                     levels = c("unexposed",
                                                "exposed (unemployed at t1)"))
pair_cc_ps$exposure2 <- factor(pair_cc_ps$exposure2,
                                     levels = c("unexposed",
                                                "exposed (job loss between t0 and t1"))


### convert SF-12 outcomes to numeric to allow svyglm to work
pair_cc_ps$sf12pcs_dv_t1 <- as.numeric(pair_cc_ps$sf12pcs_dv_t1)
pair_cc_ps$sf12mcs_dv_t1 <- as.numeric(pair_cc_ps$sf12mcs_dv_t1)


### create weighted analytic df ----------------
## add ps weights onto analytic df
#mw_spine <- pair_cc_ps %>% 
#  dplyr::select(pidp, mw)
#
#pair_cc_ps <- pair_cc_ps %>% 
#  full_join(mw_spine)
#
#pair_cc_analytic_svy <- svydesign(ids = ~1,
#                                  data = pair_cc_analytic,
#                                  weights = ~mw)
#

################################################################################
#####                           descriptive tables                         #####
################################################################################

#### unweighted outcomes table -------------------------------------------------
table_outcomes_unweighted <- CreateTableOne(vars = outcome_vector2, 
                                            strata = "exposure1",
                                            data = pair_cc_ps,
                                            test = TRUE)

table_outcomes_unweighted_sav <- print(table_outcomes_unweighted, 
                                       showAllLevels = TRUE, 
                                       smd = TRUE, 
                                       formatOptions = list(big.mark = ","))


#### weighted outcomes table ---------------------------------------------------
table_outcomes_weighted <- svyCreateTableOne(vars = outcome_vector2,
                                             strata = "exposure1",
                                             data = pair_cc_ps_svy,
                                             test = TRUE)

table_outcomes_weighted_sav <- print(table_outcomes_weighted, 
                                     showAllLevels = TRUE,  
                                     smd = TRUE)

write.csv(table_outcomes_weighted_sav, "./output/weighted_descriptives/table_outcomes_weighted_sav.csv")

################################################################################
#####                      weighted regression models                     ######
################################################################################

# initial models prior to doubly robust estimation 

### SF-12 PCS -----------------------

pcs_svyglm_mod <- svyglm(sf12pcs_dv_t1 ~ exposure1,
                         #                         family = ,
                         design = pair_cc_ps_svy, 
                         na.action = na.omit)

pcs_svyglm_summary <- summary(pcs_svyglm_mod)

pcs_svyglm_df <- tidy(pcs_svyglm_mod)

## confidence intervals
pcs_svyglm_df_ci <- data.frame(confint(pcs_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
pcs_svyglm_df_ci <- cbind(rownames(pcs_svyglm_df_ci),pcs_svyglm_df_ci, row.names=NULL)

pcs_svyglm_df_ci <- pcs_svyglm_df_ci %>% 
  rename(term = `rownames(pcs_svyglm_df_ci)`)

## join dfs together
pcs_svyglm_df <- pcs_svyglm_df %>% 
  left_join(pcs_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### SF-12 MCS -----------------------

mcs_svyglm_mod <- svyglm(sf12mcs_dv_t1 ~ exposure1,
                         #                         family = ,
                         design =pair_cc_ps_svy, 
                         na.action = na.omit)

mcs_svyglm_summary <- summary(mcs_svyglm_mod)

mcs_svyglm_df <- tidy(mcs_svyglm_mod)

## confidence intervals
mcs_svyglm_df_ci <- data.frame(confint(mcs_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) 

# add in row names
mcs_svyglm_df_ci <- cbind(rownames(mcs_svyglm_df_ci),mcs_svyglm_df_ci, row.names=NULL)

mcs_svyglm_df_ci <- mcs_svyglm_df_ci %>% 
  rename(term = `rownames(mcs_svyglm_df_ci)`)

## join dfs together
mcs_svyglm_df <- mcs_svyglm_df %>% 
  left_join(mcs_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### poor self-rated health -------------------

srh_svyglm_mod <- svyglm(srh_bin_t1 ~ exposure1,
                         family = quasibinomial,
                         design = pair_cc_ps_svy, 
                         na.action = na.omit)

srh_svyglm_summary <- summary(srh_svyglm_mod)


## coefficients dataframe
srh_svyglm_df <- tidy(srh_svyglm_mod) %>% 
  # exponentiate to get ORs
  mutate(estimate = exp(estimate))

## confidence intervals
srh_svyglm_df_ci <- data.frame(confint(srh_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
srh_svyglm_df_ci <- cbind(rownames(srh_svyglm_df_ci),srh_svyglm_df_ci, row.names=NULL)

srh_svyglm_df_ci <- srh_svyglm_df_ci %>% 
  rename(term = `rownames(srh_svyglm_df_ci)`)

## join dfs together
srh_svyglm_df <- srh_svyglm_df %>% 
  left_join(srh_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

### GHQ-12 caseness -----------------

ghq_svyglm_mod <- svyglm(ghq_case4_t1 ~ exposure1,
                         family = quasibinomial,
                         design = pair_cc_ps_svy, 
                         na.action = na.omit)

ghq_svyglm_summary <- summary(ghq_svyglm_mod)


## coefficients dataframe
ghq_svyglm_df <- tidy(ghq_svyglm_mod) %>% 
  mutate(estimate = exp(estimate))

## confidence intervals
ghq_svyglm_df_ci <- data.frame(confint(ghq_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
ghq_svyglm_df_ci <- cbind(rownames(ghq_svyglm_df_ci),ghq_svyglm_df_ci, row.names=NULL)

ghq_svyglm_df_ci <- ghq_svyglm_df_ci %>% 
  rename(term = `rownames(ghq_svyglm_df_ci)`)

## join dfs together
ghq_svyglm_df <- ghq_svyglm_df %>% 
  left_join(ghq_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))


#### combine into single dataframe

reg_df <- srh_svyglm_df %>% 
  bind_rows(ghq_svyglm_df,
            pcs_svyglm_df,
            mcs_svyglm_df) %>% 
  filter(term!="(Intercept)") %>% 
  dplyr::select(-term) %>% 
  dplyr::select(outcome, everything()) %>% 
  rename(t_value=statistic)

################################################################################
#####               double robust weighted regression models               #####
################################################################################

# remove all unneed df's
rm(ghq_svyglm_df, ghq_svyglm_df_ci, ghq_svyglm_mod, ghq_svyglm_summary,
   mcs_svyglm_df, mcs_svyglm_df_ci, mcs_svyglm_mod, mcs_svyglm_summary,
#   mw_spine,
#   pair_cc_ps,
#   pair_cc_ps_svy,
   pcs_svyglm_df, pcs_svyglm_df_ci, pcs_svyglm_mod, pcs_svyglm_summary,
#   reg_df,
   srh_svyglm_df, srh_svyglm_df_ci, srh_svyglm_mod, srh_svyglm_summary,
   table_outcomes_unweighted,table_outcomes_unweighted_sav,
   table_outcomes_weighted, table_outcomes_weighted_sav)


#### SF-12 PCS -----------------------------------------------------------------

dr_pcs_svyglm_mod <- svyglm(sf12pcs_dv_t1 ~ exposure1 +
                              sex_dv_t0 +
                              age_dv_t0 +
                              age_dv_t1 +
                              non_white_t0  +
                              marital_status_t0 +
                              marital_status_t1 +
                              hiqual_dv_t0 +
                              gor_dv_t0 +
                              gor_dv_t1 +
                              sic2007_section_lab_t0 +
#                              sic2007_section_lab_t1 +
                              soc2000_major_group_title_t0 +
#                              soc2000_major_group_title_t1 +
                              jbft_dv_t0 +
#                              jbft_dv_t1 +
                              small_firm_t0 +
#                              small_firm_t1 +
#                             jbhrs_t0 +
#                              jbhrs_t1 +
                              emp_contract_t0 +
#                              emp_contract_t1 +
                              broken_emp_t0 +
                              j2has_dv_t0 +
#                              j2has_dv_t1 +
                              rel_pov_t0 +
                              health_t0 +
                              health_t1 +
                              sf12pcs_dv_t0,
#                           family = ,
                            design = pair_cc_ps_svy, 
                            na.action = na.omit)

dr_pcs_svyglm_summary <- summary(dr_pcs_svyglm_mod)

dr_pcs_svyglm_df <- tidy(dr_pcs_svyglm_mod)

## confidence intervals
dr_pcs_svyglm_df_ci <- data.frame(confint(dr_pcs_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)

# add in row names
dr_pcs_svyglm_df_ci <- cbind(rownames(dr_pcs_svyglm_df_ci),dr_pcs_svyglm_df_ci, row.names=NULL)

dr_pcs_svyglm_df_ci <- dr_pcs_svyglm_df_ci %>% 
  rename(term = `rownames(dr_pcs_svyglm_df_ci)`)

## join dfs together
dr_pcs_svyglm_df <- dr_pcs_svyglm_df %>% 
  left_join(dr_pcs_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 PCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### SF-12 MCS -----------------------------------------------------------------

dr_mcs_svyglm_mod <- svyglm(sf12mcs_dv_t1 ~ exposure1 +
                              sex_dv_t0 +
                              age_dv_t0 +
                              age_dv_t1 +
                              non_white_t0  +
                              marital_status_t0 +
                              marital_status_t1 +
                              hiqual_dv_t0 +
                              gor_dv_t0 +
                              gor_dv_t1 +
                              sic2007_section_lab_t0 +
#                              sic2007_section_lab_t1 +
                              soc2000_major_group_title_t0 +
#                              soc2000_major_group_title_t1 +
                              jbft_dv_t0 +
#                              jbft_dv_t1 +
                              small_firm_t0 +
#                              small_firm_t1 +
#                              jbhrs_t0 +
#                              jbhrs_t1 +
                              emp_contract_t0 +
#                              emp_contract_t1 +
                              broken_emp_t0 +
                              j2has_dv_t0 +
#                              j2has_dv_t1 +
                              rel_pov_t0 +
                              health_t0 +
                              health_t1 +
                              sf12mcs_dv_t0,
#                         family = ,
                            design = pair_cc_ps_svy, 
                            na.action = na.omit)

dr_mcs_svyglm_summary <- summary(dr_mcs_svyglm_mod)

dr_mcs_svyglm_df <- tidy(dr_mcs_svyglm_mod)

## confidence intervals
dr_mcs_svyglm_df_ci <- data.frame(confint(dr_mcs_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) #%>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))

# add in row names
dr_mcs_svyglm_df_ci <- cbind(rownames(dr_mcs_svyglm_df_ci),dr_mcs_svyglm_df_ci, row.names=NULL)

dr_mcs_svyglm_df_ci <- dr_mcs_svyglm_df_ci %>% 
  rename(term = `rownames(dr_mcs_svyglm_df_ci)`)

## join dfs together
dr_mcs_svyglm_df <- dr_mcs_svyglm_df %>% 
  left_join(dr_mcs_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "SF-12 MCS",
         est_type = "coefficient",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### poor self-rated health ----------------------------------------------------

dr_srh_svyglm_mod <- svyglm(srh_bin_t1 ~ exposure1 +
                              sex_dv_t0 +
                              age_dv_t0 +
                              age_dv_t1 +
                              non_white_t0  +
                              marital_status_t0 +
                              marital_status_t1 +
                              hiqual_dv_t0 +
                              gor_dv_t0 +
                              gor_dv_t1 +
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
                              srh_bin_t0,
                            family = quasibinomial,
                            design = pair_cc_ps_svy, 
                            na.action = na.omit)
## note - model not currently converging due to coding error on srh_bin_t0 == re-run once PS code done

dr_srh_svyglm_summary <- summary(dr_srh_svyglm_mod)


## coefficients dataframe
dr_srh_svyglm_df <- tidy(dr_srh_svyglm_mod) %>% 
  # exponentiate to get ORs
  mutate(estimate = exp(estimate))

## confidence intervals
dr_srh_svyglm_df_ci <- data.frame(confint(dr_srh_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_srh_svyglm_df_ci <- cbind(rownames(dr_srh_svyglm_df_ci),dr_srh_svyglm_df_ci, row.names=NULL)

dr_srh_svyglm_df_ci <- dr_srh_svyglm_df_ci %>% 
  rename(term = `rownames(dr_srh_svyglm_df_ci)`)

## join dfs together
dr_srh_svyglm_df <- dr_srh_svyglm_df %>% 
  left_join(dr_srh_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "Poor self-rated health",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))

#### GHQ-12 caseness -----------------------------------------------------------

dr_ghq_svyglm_mod <- svyglm(ghq_case4_t1 ~ exposure1 +
                              sex_dv_t0 +
                              age_dv_t0 +
                              age_dv_t1 +
                              non_white_t0  +
                              marital_status_t0 +
                              marital_status_t1 +
                              hiqual_dv_t0 +
                              gor_dv_t0 +
                              gor_dv_t1 +
                              sic2007_section_lab_t0 +
#                              sic2007_section_lab_t1 +
                              soc2000_major_group_title_t0 +
#                              soc2000_major_group_title_t1 +
                              jbft_dv_t0 +
#                              jbft_dv_t1 +
                              small_firm_t0 +
#                              small_firm_t1 +
#                              jbhrs_t0 +
#                              jbhrs_t1 +
                              emp_contract_t0 +
#                              emp_contract_t1 +
                              broken_emp_t0 +
                              j2has_dv_t0 +
#                              j2has_dv_t1 +
                              rel_pov_t0 +
                              health_t0 +
                              health_t1 +
                              ghq_case4_t0,
                            family = quasibinomial,
                            design = pair_cc_ps_svy, 
                            na.action = na.omit)

dr_ghq_svyglm_summary <- summary(dr_ghq_svyglm_mod)


## coefficients dataframe
dr_ghq_svyglm_df <- tidy(dr_ghq_svyglm_mod) %>% 
  mutate(estimate = exp(estimate))

## confidence intervals
dr_ghq_svyglm_df_ci <- data.frame(confint(dr_ghq_svyglm_mod)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..) %>% 
  mutate(lci = exp(lci),
         uci = exp(uci))

# add in row names
dr_ghq_svyglm_df_ci <- cbind(rownames(dr_ghq_svyglm_df_ci),dr_ghq_svyglm_df_ci, row.names=NULL)

dr_ghq_svyglm_df_ci <- dr_ghq_svyglm_df_ci %>% 
  rename(term = `rownames(dr_ghq_svyglm_df_ci)`)

## join dfs together
dr_ghq_svyglm_df <- dr_ghq_svyglm_df %>% 
  left_join(dr_ghq_svyglm_df_ci) %>% 
  mutate(term = str_remove(term, "exposure1"),
         outcome = "GHQ-12 caseness (4+)",
         est_type = "OR",
         p.value = ifelse(p.value<0.001,"<0.001",
                          ifelse(p.value<0.01,"<0.01",
                                 ifelse(p.value<0.05,"<0.05",       
                                        p.value))))


#### combine into single dataframe

dr_reg_df <- dr_pcs_svyglm_df %>% 
  bind_rows(dr_mcs_svyglm_df,
            dr_srh_svyglm_df,
            dr_ghq_svyglm_df) %>% 
  filter(term=="exposed (unemployed at t1)") %>% 
  dplyr::select(-term) %>% 
  dplyr::select(outcome, everything()) %>% 
  rename(t_value=statistic)

write.csv(table_outcomes_weighted_sav, "./output/weighted_outcomes/cc_double_robust_MSM.csv")
