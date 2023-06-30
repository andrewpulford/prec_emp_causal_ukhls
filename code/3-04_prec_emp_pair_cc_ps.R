################################################################################

# Precarious employment and health - Understanding Society
# 3-04 - Paired unweighted complete case analysis propensity scores for risk of 
# job loss 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) produces unweighted risk of job loss propensity scores for complete case 
#     paired data


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
library(lme4) # for multi-level modelling

################################################################################
#####                         load and prepare data                        #####
################################################################################

####load eligible cases --------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")

### convert binary outcome vars to factors to allow svyglm to work
pair_cc_analytic$srh_bin_t1 <- factor(pair_cc_analytic$srh_bin_t0)
pair_cc_analytic$ghq_case4_t1 <- factor(pair_cc_analytic$ghq_case4_t1)

### convert SF-12 outcomes to numeric to allow svyglm to work
pair_cc_analytic$sf12pcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t1)
pair_cc_analytic$sf12mcs_dv_t1 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t1)


#### create vectors for variables ----------------------------------------------
id_wt_vector <- c("pidp", "psu", "strata", "wt_name", "wt_value")

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
                "fimnnet_dv_t0",
                "jbhrs_t0",
                "health_t0",
                "srh_bin_t0",
                "ghq_case4_t0",
                "sf12mcs_dv_t0",
                "sf12pcs_dv_t0")
# use missing cat
# don't create separate outcome df's unless high # missing

cov_vector2 <- c("age_dv_t1",  
                "marital_status_t1",
                "gor_dv_t1",
                "sic2007_section_lab_t1",
                "soc2000_major_group_title_t1",
                "jbft_dv_t1",
                "small_firm_t1",
                "emp_contract_t1",
                "broken_emp_t1",
                "j2has_dv_t1",
                "fimnnet_dv_t1",
                "jbhrs_t1",
                "health_t1",
                "srh_bin_t1",
                "ghq_case4_t1",
                "sf12mcs_dv_t1",
                "sf12pcs_dv_t1")


outcome_vector <- c("srh_bin_t1",
                    "ghq_case4_t1",
                    "sf12mcs_dv_t1",
                    "sf12pcs_dv_t1",
                    "exposure1",
                    "exposure2")

outcome_vector2 <- c("srh_bin_t1",
                    "ghq_case4_t1",
                    "sf12mcs_dv_t1",
                    "sf12pcs_dv_t1")

#### keep only variables required for propensity score -------------------------
pair_cc_ps <- pair_cc_analytic %>% 
  dplyr::select(all_of(c(id_wt_vector, cov_vector, outcome_vector))) %>% 
  dplyr::select(-c(psu, strata, wt_name, wt_value))


#### convert scale vars to numeric ---------------------------------------------
pair_cc_ps <- pair_cc_ps %>% 
  mutate(fimnnet_dv_t0 = as.numeric(fimnnet_dv_t0)) %>%
  ## these to character first so score converted rather than factor level
  mutate(sf12mcs_dv_t0 = as.character(sf12mcs_dv_t0),
         sf12pcs_dv_t0 = as.character(sf12pcs_dv_t0),
         sf12mcs_dv_t1 = as.character(sf12mcs_dv_t1),
         sf12pcs_dv_t1 = as.character(sf12pcs_dv_t1)) %>% 
  mutate(sf12mcs_dv_t0 = as.numeric(sf12mcs_dv_t0),
         sf12pcs_dv_t0 = as.numeric(sf12pcs_dv_t0),
         sf12mcs_dv_t1 = as.numeric(sf12mcs_dv_t1),
         sf12pcs_dv_t1 = as.numeric(sf12pcs_dv_t1))


pair_cc_ps[sapply(pair_cc_ps, is.character)] <- lapply(pair_cc_ps[sapply(pair_cc_ps, is.character)], as.factor)

#### revel exposure vars so that unexposed is lower ref category ---------------
pair_cc_ps$exposure1 <- factor(pair_cc_ps$exposure1,
                                     levels = rev(levels(pair_cc_ps$exposure1)))

pair_cc_ps$exposure2 <- factor(pair_cc_ps$exposure2,
                                     levels = rev(levels(pair_cc_ps$exposure2)))


#### remove any NAs not coded as a missing category ----------------------------
##  check NAs
sapply(pair_cc_ps, function(x) sum(is.na(x)))
## remove
pair_cc_ps <- na.omit(pair_cc_ps)
## check again
sapply(pair_cc_ps, function(x) sum(is.na(x)))


### use this if full df is too big
#test_df <- pair_cc_ps %>% 
#  slice_sample(prop = 1)

#### load unmatched SMD df
table_one_unmatched_smd <- read.csv("./working_data/table_one_unmatched_smd.csv") %>% 
  dplyr::select(c(var, smd, imbalance_flag, matched))


################################################################################
#####                               functions                              #####
################################################################################

#### propensity score model ----------------------------------------------------

ps_model <- function(data = pair_cc_ps, outcome){
  glmer(outcome ~
         sex_dv_t0 +
         age_dv_t0 +
         non_white_t0  +
         marital_status_t0 +
         hiqual_dv_t0 +
         gor_dv_t0 +
         sic2007_section_lab_t0 +
         soc2000_major_group_title_t0 +
         jbhrs_t0 +
         emp_contract_t0 +
         broken_emp_t0 +
         j2has_dv_t0 +
         fimnnet_dv_t0 +
         health_t0 +
         srh_bin_t0 +
         ghq_case4_t0 +
         sf12mcs_dv_t0 +
         sf12pcs_dv_t0 +
         1+(1|pidp),
            family = binomial(link="logit"),
            data = data)
  
  
}

#### outcome model -------------------------------------------------------------

outcome_model <- function(outcome, exposure = exposure1, data = pair_cc_ps_svy){
  svyglm(outcome ~ exposure,
      family = quasibinomial,
      design = data, 
      na.action = na.omit)
  
  
}


################################################################################
#####                     Unemployment at t1 PS model                      #####
################################################################################

### call the function
ps_mod_exp1 <- ps_model(data = pair_cc_ps, outcome = pair_cc_ps$exposure1)

### summary of model
summary(ps_mod_exp1)


### predicted probability of being assigned to exposed group
pair_cc_ps$ps_exp1 <- predict(ps_mod_exp1, type = "response")
summary(pair_cc_ps$ps_exp1)

### predicted probability of being assigned to unexposed group
pair_cc_ps$ps_noexp1 <- 1-pair_cc_ps$ps_exp1
summary(pair_cc_ps$ps_noexp1)

### predicted probability of being assigned to actual exposure status
pair_cc_ps$ps_assign <- NA
pair_cc_ps$ps_assign[pair_cc_ps$exposure1=="exposed (unemployed at t1)"] <- pair_cc_ps$ps_exp1[pair_cc_ps$exposure1=="exposed (unemployed at t1)"]
pair_cc_ps$ps_assign[pair_cc_ps$exposure1=="unexposed"] <- pair_cc_ps$ps_noexp1[pair_cc_ps$exposure1=="unexposed"]

### smaller of ps_exp1 and ps_noexp1 for matchnig weight

pair_cc_ps$ps_min <- pmin(pair_cc_ps$ps_exp1, pair_cc_ps$ps_noexp1)

################################################################################
#####                       propensity score matching                      #####
################################################################################

### remove missing as category for tables == don't use for now, throws error
#pair_cc_ps <- pair_cc_ps %>% 
#  mutate(across(.cols = everything(), 
#                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 


#### match propensity scores using Matching package ----------------------------
list_match <- Match(Tr = (pair_cc_ps$exposure1=="exposed (unemployed at t1)"),
# logit of PS/1-PS
X = log(pair_cc_ps$ps_exp1/pair_cc_ps$ps_noexp1),
## 1:1 matching ratio
M = 1,
## caliper = 0.2 * SD(logit(PS))
caliper  = 0.2,
replace  = FALSE,
ties     = TRUE,
version  = "fast")

#### extract matched data ------------------------------------------------------
df_matched <- pair_cc_ps[unlist(list_match[c("index.treated","index.control")]), ]

#### matched Table One with SMD ------------------------------------------------
table_one_matched <- CreateTableOne(vars = cov_vector,
                                    strata = "exposure1",
                                    data = df_matched,
                                    test = FALSE)

table_one_matched_sav <- print(table_one_matched,  smd = TRUE)

write.csv(table_one_matched_sav, "./output/matched_descriptives/table_one_matched.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(table_one_matched) > 0.1))


################################################################################
#####                   propensity score matching weight                   #####
################################################################################

### create matching weight
pair_cc_ps$mw <- pair_cc_ps$ps_min/pair_cc_ps$ps_assign

### create weighted data
pair_cc_ps_svy <- svydesign(ids = ~1,
                                  data = pair_cc_ps,
                                  weights = ~mw)

### weighted table one
table_one_weighted <- svyCreateTableOne(vars = cov_vector,
                                        strata = "exposure1",
                                        data = pair_cc_ps_svy,
                                        test = FALSE)

table_one_weighted_sav <- print(table_one_weighted,  smd = TRUE)

write.csv(table_one_weighted_sav, "./output/weighted_descriptives/table_one_weighted.csv")

### count covariates with an important imbalance (>0.1)



########## propensity score overlap weight??????????????????????????????????????


#### write working files -------------------------------------------------------
write_rds(pair_cc_ps, "working_data/pair_cc_ps.rds")


################################################################################
#####                  assess balance pre/post-matching                   ######
################################################################################

#### create df for matched smd values ------------------------------------------
table_one_matched_smd <- data.frame(ExtractSmd(table_one_matched))
table_one_matched_smd <- table_one_matched_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "matched")

#### create df for weighted smd values -----------------------------------------
table_one_weighted_smd <- data.frame(ExtractSmd(table_one_weighted))
table_one_weighted_smd <- table_one_weighted_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "weighted")

#### bind to unmatched smd df --------------------------------------------------
assess_matching_balance <- table_one_matched_smd %>% 
  bind_rows(table_one_weighted_smd,table_one_unmatched_smd)

### plot
unemp_t1_balance_plot <- assess_matching_balance %>% 
  ggplot(aes(x=smd, y=var, col=matched)) +
  geom_jitter(height = 0.01, width = 0.01) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("blue","red", "green"))

tiff("./output/weighted_descriptives/unemp_t1_balance_plot.tiff")
unemp_t1_balance_plot
dev.off()

#################################################################################
######                           outcome analysis                           #####
#################################################################################
#
##### unweighted outcomes table -------------------------------------------------
#table_outcomes_unweighted <- CreateTableOne(vars = outcome_vector2, 
#                                            strata = "exposure1",
#                                            data = pair_cc_ps,
#                                            test = TRUE)
#
#table_outcomes_unweighted_sav <- print(table_outcomes_unweighted, 
#                                       showAllLevels = TRUE, 
#                                       smd = TRUE, 
#                       formatOptions = list(big.mark = ","))
#
#
#### weighted outcomes table
#table_outcomes_weighted <- svyCreateTableOne(vars = outcome_vector2,
#                                        strata = "exposure1",
#                                        data = pair_cc_ps_svy,
#                                        test = TRUE)
#
#table_outcomes_weighted_sav <- print(table_outcomes_weighted, 
#                                     showAllLevels = TRUE,  
#                                     smd = TRUE)
#
#write.csv(table_outcomes_weighted_sav, "./output/weighted_descriptives/table_outcomes_weighted_sav.csv")
#
##### weighted regression models ------------------------------------------------
#
#### create weighted analytic df ----------------
### add ps weights onto analytic df
#mw_spine <- pair_cc_ps %>% 
#  dplyr::select(pidp, mw)
#
#pair_cc_analytic <- pair_cc_analytic %>% 
#  right_join(mw_spine)
#
#pair_cc_analytic_svy <- svydesign(ids = ~1,
#                            data = pair_cc_analytic,
#                            weights = ~mw)
#
#
#### poor self-rated health -------------------
#
#srh_svyglm_mod <- svyglm(srh_bin_t1 ~ exposure1,
#       family = quasibinomial,
#       design = pair_cc_analytic_svy, 
#       na.action = na.omit)
#
#srh_svyglm_summary <- summary(srh_svyglm_mod)
#
#
### coefficients dataframe
#srh_svyglm_df <- tidy(srh_svyglm_mod) %>% 
## exponentiate to get ORs
#  mutate(estimate = exp(estimate))
#
### confidence intervals
#srh_svyglm_df_ci <- data.frame(confint(srh_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#srh_svyglm_df_ci <- cbind(rownames(srh_svyglm_df_ci),srh_svyglm_df_ci, row.names=NULL)
#
#srh_svyglm_df_ci <- srh_svyglm_df_ci %>% 
#  rename(term = `rownames(srh_svyglm_df_ci)`)
#
### join dfs together
#srh_svyglm_df <- srh_svyglm_df %>% 
#  left_join(srh_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "Poor self-rated health",
#         est_type = "OR",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
#### GHQ-12 caseness -----------------
#
#ghq_svyglm_mod <- svyglm(ghq_case4_t1 ~ exposure1,
#                         family = quasibinomial,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#ghq_svyglm_summary <- summary(ghq_svyglm_mod)
#
#
### coefficients dataframe
#ghq_svyglm_df <- tidy(ghq_svyglm_mod) %>% 
#  mutate(estimate = exp(estimate))
#
### confidence intervals
#ghq_svyglm_df_ci <- data.frame(confint(ghq_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#ghq_svyglm_df_ci <- cbind(rownames(ghq_svyglm_df_ci),ghq_svyglm_df_ci, row.names=NULL)
#
#ghq_svyglm_df_ci <- ghq_svyglm_df_ci %>% 
#  rename(term = `rownames(ghq_svyglm_df_ci)`)
#
### join dfs together
#ghq_svyglm_df <- ghq_svyglm_df %>% 
#  left_join(ghq_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "GHQ-12 caseness (4+)",
#         est_type = "OR",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
#### SF-12 PCS -----------------------
#
#pcs_svyglm_mod <- svyglm(sf12pcs_dv_t1 ~ exposure1,
##                         family = ,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#pcs_svyglm_summary <- summary(pcs_svyglm_mod)
#
#pcs_svyglm_df <- tidy(pcs_svyglm_mod)
#
### confidence intervals
#pcs_svyglm_df_ci <- data.frame(confint(pcs_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#pcs_svyglm_df_ci <- cbind(rownames(pcs_svyglm_df_ci),pcs_svyglm_df_ci, row.names=NULL)
#
#pcs_svyglm_df_ci <- pcs_svyglm_df_ci %>% 
#  rename(term = `rownames(pcs_svyglm_df_ci)`)
#
### join dfs together
#pcs_svyglm_df <- pcs_svyglm_df %>% 
#  left_join(pcs_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "SF-12 PCS",
#         est_type = "coefficient",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
#### SF-12 MCS -----------------------
#
#mcs_svyglm_mod <- svyglm(sf12mcs_dv_t1 ~ exposure1,
#                         #                         family = ,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#mcs_svyglm_summary <- summary(mcs_svyglm_mod)
#
#mcs_svyglm_df <- tidy(mcs_svyglm_mod)
#
### confidence intervals
#mcs_svyglm_df_ci <- data.frame(confint(mcs_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#mcs_svyglm_df_ci <- cbind(rownames(mcs_svyglm_df_ci),mcs_svyglm_df_ci, row.names=NULL)
#
#mcs_svyglm_df_ci <- mcs_svyglm_df_ci %>% 
#  rename(term = `rownames(mcs_svyglm_df_ci)`)
#
### join dfs together
#mcs_svyglm_df <- mcs_svyglm_df %>% 
#  left_join(mcs_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "SF-12 MCS",
#         est_type = "coefficient",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                          p.value))))
#
##### combine into single dataframe
#
#reg_df <- srh_svyglm_df %>% 
#  bind_rows(ghq_svyglm_df,
#            pcs_svyglm_df,
#            mcs_svyglm_df) %>% 
#  filter(term!="(Intercept)") %>% 
#  dplyr::select(-term) %>% 
#  dplyr::select(outcome, everything()) %>% 
#  rename(t_value=statistic)
#
##### doubly-robust weighted regression models ----------------------------------
#
#### poor self-rated health -------------------
#
#dr_srh_svyglm_mod <- svyglm(srh_bin_t1 ~ exposure1 +
#                           sex_dv_t0 +
#                           age_dv_t0 +
#                             age_dv_t1 +
#                           non_white_t0  +
#                           marital_status_t0 +
#                             marital_status_t1 +
#                           hiqual_dv_t0 +
#                             hiqual_dv_t1 +
#                           gor_dv_t0 +
#                             gor_dv_t1 +
#                           sic2007_section_lab_t0 +
#                             sic2007_section_lab_t1 +
#                           soc2000_major_group_title_t0 +
#                             soc2000_major_group_title_t1 +
#                           jbhrs_t0 +
#                             jbhrs_t1 +
#                           emp_contract_t0 +
#                             emp_contract_t1 +
#                           broken_emp_t0 +
#                             broken_emp_t1 +
#                           j2has_dv_t0 +
#                             j2has_dv_t1 +
#                           fimnnet_dv_t0 +
#                             fimnnet_dv_t1 +
#                           health_t0 +
#                             health_t1 +
#                           srh_bin_t0 +
#                           ghq_case4_t0 +
#                             ghq_case4_t1 +
#                           sf12mcs_dv_t0 +
#                             sf12mcs_dv_t1 +
#                             sf12pcs_dv_t0 +
#                           sf12pcs_dv_t1,
#                           family = quasibinomial,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#dr_srh_svyglm_summary <- summary(dr_srh_svyglm_mod)
#
#
### coefficients dataframe
#dr_srh_svyglm_df <- tidy(dr_srh_svyglm_mod) %>% 
#  # exponentiate to get ORs
#  mutate(estimate = exp(estimate))
#
### confidence intervals
#dr_srh_svyglm_df_ci <- data.frame(confint(dr_srh_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#dr_srh_svyglm_df_ci <- cbind(rownames(dr_srh_svyglm_df_ci),dr_srh_svyglm_df_ci, row.names=NULL)
#
#dr_srh_svyglm_df_ci <- dr_srh_svyglm_df_ci %>% 
#  rename(term = `rownames(dr_srh_svyglm_df_ci)`)
#
### join dfs together
#dr_srh_svyglm_df <- dr_srh_svyglm_df %>% 
#  left_join(dr_srh_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "Poor self-rated health",
#         est_type = "OR",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
#### GHQ-12 caseness -----------------
#
#dr_ghq_svyglm_mod <- svyglm(ghq_case4_t1 ~ exposure1 +
#                         sex_dv_t0 +
#                           age_dv_t0 +
#                           non_white_t0  +
#                           marital_status_t0 +
#                           hiqual_dv_t0 +
#                           gor_dv_t0 +
#                           sic2007_section_lab_t0 +
#                           soc2000_major_group_title_t0 +
#                           jbhrs_t0 +
#                           emp_contract_t0 +
#                           broken_emp_t0 +
#                           j2has_dv_t0 +
#                           fimnnet_dv_t0 +
#                           health_t0 +
#                           srh_bin_t0 +
#                           ghq_case4_t0 +
#                           sf12mcs_dv_t0 +
#                           sf12pcs_dv_t0,
#                         family = quasibinomial,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#dr_ghq_svyglm_summary <- summary(dr_ghq_svyglm_mod)
#
#
### coefficients dataframe
#dr_ghq_svyglm_df <- tidy(dr_ghq_svyglm_mod) %>% 
#  mutate(estimate = exp(estimate))
#
### confidence intervals
#dr_ghq_svyglm_df_ci <- data.frame(confint(dr_ghq_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#dr_ghq_svyglm_df_ci <- cbind(rownames(dr_ghq_svyglm_df_ci),dr_ghq_svyglm_df_ci, row.names=NULL)
#
#dr_ghq_svyglm_df_ci <- dr_ghq_svyglm_df_ci %>% 
#  rename(term = `rownames(dr_ghq_svyglm_df_ci)`)
#
### join dfs together
#dr_ghq_svyglm_df <- dr_ghq_svyglm_df %>% 
#  left_join(dr_ghq_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "GHQ-12 caseness (4+)",
#         est_type = "OR",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
#### SF-12 PCS -----------------------
#
#dr_pcs_svyglm_mod <- svyglm(sf12pcs_dv_t1 ~ exposure1 +
#                         sex_dv_t0 +
#                           age_dv_t0 +
#                           non_white_t0  +
#                           marital_status_t0 +
#                           hiqual_dv_t0 +
#                           gor_dv_t0 +
#                           sic2007_section_lab_t0 +
#                           soc2000_major_group_title_t0 +
#                           jbhrs_t0 +
#                           emp_contract_t0 +
#                           broken_emp_t0 +
#                           j2has_dv_t0 +
#                           fimnnet_dv_t0 +
#                           health_t0 +
#                           srh_bin_t0 +
#                           ghq_case4_t0 +
#                           sf12mcs_dv_t0 +
#                           sf12pcs_dv_t0,
#                         #                         family = ,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#dr_pcs_svyglm_summary <- summary(dr_pcs_svyglm_mod)
#
#dr_pcs_svyglm_df <- tidy(dr_pcs_svyglm_mod)
#
### confidence intervals
#dr_pcs_svyglm_df_ci <- data.frame(confint(dr_pcs_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#dr_pcs_svyglm_df_ci <- cbind(rownames(dr_pcs_svyglm_df_ci),dr_pcs_svyglm_df_ci, row.names=NULL)
#
#dr_pcs_svyglm_df_ci <- dr_pcs_svyglm_df_ci %>% 
#  rename(term = `rownames(dr_pcs_svyglm_df_ci)`)
#
### join dfs together
#dr_pcs_svyglm_df <- dr_pcs_svyglm_df %>% 
#  left_join(dr_pcs_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "SF-12 PCS",
#         est_type = "coefficient",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
#### SF-12 MCS -----------------------
#
#dr_mcs_svyglm_mod <- svyglm(sf12mcs_dv_t1 ~ exposure1 +
#                         sex_dv_t0 +
#                           age_dv_t0 +
#                           non_white_t0  +
#                           marital_status_t0 +
#                           hiqual_dv_t0 +
#                           gor_dv_t0 +
#                           sic2007_section_lab_t0 +
#                           soc2000_major_group_title_t0 +
#                           jbhrs_t0 +
#                           emp_contract_t0 +
#                           broken_emp_t0 +
#                           j2has_dv_t0 +
#                           fimnnet_dv_t0 +
#                           health_t0 +
#                           srh_bin_t0 +
#                           ghq_case4_t0 +
#                           sf12mcs_dv_t0 +
#                           sf12pcs_dv_t0,
#                         #                         family = ,
#                         design = pair_cc_analytic_svy, 
#                         na.action = na.omit)
#
#dr_mcs_svyglm_summary <- summary(dr_mcs_svyglm_mod)
#
#dr_mcs_svyglm_df <- tidy(dr_mcs_svyglm_mod)
#
### confidence intervals
#dr_mcs_svyglm_df_ci <- data.frame(confint(dr_mcs_svyglm_mod)) %>% 
#  rename(lci = X2.5..,
#         uci = X97.5..) %>% 
#  mutate(lci = exp(lci),
#         uci = exp(uci))
#
## add in row names
#dr_mcs_svyglm_df_ci <- cbind(rownames(dr_mcs_svyglm_df_ci),dr_mcs_svyglm_df_ci, row.names=NULL)
#
#dr_mcs_svyglm_df_ci <- dr_mcs_svyglm_df_ci %>% 
#  rename(term = `rownames(dr_mcs_svyglm_df_ci)`)
#
### join dfs together
#dr_mcs_svyglm_df <- dr_mcs_svyglm_df %>% 
#  left_join(dr_mcs_svyglm_df_ci) %>% 
#  mutate(term = str_remove(term, "exposure1"),
#         outcome = "SF-12 MCS",
#         est_type = "coefficient",
#         p.value = ifelse(p.value<0.001,"<0.001",
#                          ifelse(p.value<0.01,"<0.01",
#                                 ifelse(p.value<0.05,"<0.05",       
#                                        p.value))))
#
##### combine into single dataframe
#
#dr_reg_df <- dr_srh_svyglm_df %>% 
#  bind_rows(dr_ghq_svyglm_df,
#            dr_pcs_svyglm_df,
#            dr_mcs_svyglm_df) %>% 
#  filter(term=="exposed (unemployed at t1)") %>% 
#  dplyr::select(-term) %>% 
#  dplyr::select(outcome, everything()) %>% 
#  rename(t_value=statistic)
#
#