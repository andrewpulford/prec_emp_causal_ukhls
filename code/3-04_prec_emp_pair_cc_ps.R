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


################################################################################
#####                         load and prepare data                        #####
################################################################################

####load eligible cases --------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")


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
                "emp_contract_t0",
#                "broken_emp_t0",
                "j2has_dv_t0",
                "fimnnet_dv_t0",
                "health_t0",
                "srh_bin_t0",
                "ghq_case4_t0",
                "sf12mcs_dv_t0",
                "sf12pcs_dv_t0")
# use missing cat
# don't create separate outcome df's unless high # missing

outcome_vector <- c("srh_bin_t1",
                    "ghq_case4_t1",
                    "sf12mcs_dv_t1",
                    "sf12pcs_dv_t1",
                    "exposure1",
                    "exposure2")

#### keep only variables required for analysis ---------------------------------
pair_cc_analytic <- pair_cc_analytic %>% 
  dplyr::select(c(id_wt_vector, cov_vector, outcome_vector)) %>% 
  dplyr::select(-c(psu, strata, wt_name, wt_value))


#### convert scale vars to numeric ---------------------------------------------
pair_cc_analytic <- pair_cc_analytic %>% 
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


pair_cc_analytic[sapply(pair_cc_analytic, is.character)] <- lapply(pair_cc_analytic[sapply(pair_cc_analytic, is.character)], as.factor)

#### revel exposure vars so that unexposed is lower ref category ---------------
pair_cc_analytic$exposure1 <- factor(pair_cc_analytic$exposure1,
                                     levels = rev(levels(pair_cc_analytic$exposure1)))

pair_cc_analytic$exposure2 <- factor(pair_cc_analytic$exposure2,
                                     levels = rev(levels(pair_cc_analytic$exposure2)))


#### remove any NAs not coded as a missing category ----------------------------
##  check NAs
sapply(pair_cc_analytic, function(x) sum(is.na(x)))
## remove
pair_cc_analytic <- na.omit(pair_cc_analytic)
## check again
sapply(pair_cc_analytic, function(x) sum(is.na(x)))


### use this if full df is too big
#test_df <- pair_cc_analytic %>% 
#  slice_sample(prop = 1)

#### load unmatched SMD df
table_one_unmatched_smd <- read.csv("./output/table_one_unmatched_smd.csv") %>% 
  dplyr::select(c(var, smd, imbalance_flag, matched))


################################################################################
#####                               functions                              #####
################################################################################

#### propensity score model ----------------------------------------------------

ps_model <- function(data = pair_cc_analytic, outcome){
  glm(formula = outcome ~
              sex_dv_t0 +
              age_dv_t0 +
              non_white_t0  +
              marital_status_t0 +
              hiqual_dv_t0 +
              gor_dv_t0 +
              sic2007_section_lab_t0 +
              soc2000_major_group_title_t0 +
              emp_contract_t0 +
#              broken_emp_t0 +
              j2has_dv_t0 +
              fimnnet_dv_t0 +
              health_t0 +
              srh_bin_t0 +
              ghq_case4_t0 +
              sf12mcs_dv_t0 +
              sf12pcs_dv_t0,
            family = binomial(link="logit"),
            data = data)
  
  
}


################################################################################
#####                     Unemployment at t1 PS model                      #####
################################################################################

### call the function
ps_mod_exp1 <- ps_model(data = pair_cc_analytic, outcome = pair_cc_analytic$exposure1)

### summary of model
summary(ps_mod_exp1)


### predicted probability of being assigned to exposed group
pair_cc_analytic$ps_exp1 <- predict(ps_mod_exp1, type = "response")
summary(pair_cc_analytic$ps_exp1)

### predicted probability of being assigned to unexposed group
pair_cc_analytic$ps_noexp1 <- 1-pair_cc_analytic$ps_exp1
summary(pair_cc_analytic$ps_noexp1)

### predicted probability of being assigned to actual exposure status
pair_cc_analytic$ps_assign <- NA
pair_cc_analytic$ps_assign[pair_cc_analytic$exposure1=="exposed (unemployed at t1)"] <- pair_cc_analytic$ps_exp1[pair_cc_analytic$exposure1=="exposed (unemployed at t1)"]
pair_cc_analytic$ps_assign[pair_cc_analytic$exposure1=="unexposed"] <- pair_cc_analytic$ps_noexp1[pair_cc_analytic$exposure1=="unexposed"]

### smaller of ps_exp1 and ps_noexp1 for matchnig weight

pair_cc_analytic$ps_min <- pmin(pair_cc_analytic$ps_exp1, pair_cc_analytic$ps_noexp1)

################################################################################
#####                       propensity score matching                      #####
################################################################################


#### match propensity scores using Matching package ----------------------------
list_match <- Match(Tr = (pair_cc_analytic$exposure1=="exposed (unemployed at t1)"),
# logit of PS/1-PS
X = log(pair_cc_analytic$ps_exp1/pair_cc_analytic$ps_noexp1),
## 1:1 matching ratio
M = 1,
## caliper = 0.2 * SD(logit(PS))
caliper  = 0.2,
replace  = FALSE,
ties     = TRUE,
version  = "fast")

#### extract matched data ------------------------------------------------------
df_matched <- pair_cc_analytic[unlist(list_match[c("index.treated","index.control")]), ]

#### matched Table One with SMD ------------------------------------------------
table_one_matched <- CreateTableOne(vars = cov_vector,
                                    strata = "exposure1",
                                    data = df_matched,
                                    test = FALSE)

table_one_matched_sav <- print(table_one_matched,  smd = TRUE)

write.csv(table_one_matched_sav, "./output/table_one_matched.csv")


### count covariates with an important imbalance (>0.1)
addmargins(table(ExtractSmd(table_one_matched) > 0.1))


################################################################################
#####                   propensity score matching weight                   #####
################################################################################

### create matching weight
pair_cc_analytic$mw <- pair_cc_analytic$ps_min/pair_cc_analytic$ps_assign

### create weighted data
pair_cc_analytic_svy <- svydesign(ids = ~1,
                                  data = pair_cc_analytic,
                                  weights = ~mw)

### weighted table one
table_one_weighted <- svyCreateTableOne(vars = cov_vector,
                                        strata = "exposure1",
                                        data = pair_cc_analytic_svy,
                                        test = FALSE)

table_one_weighted_sav <- print(table_one_weighted,  smd = TRUE)

write.csv(table_one_weighted_sav, "./output/table_one_weighted_sav.csv")

### count covariates with an important imbalance (>0.1)



########## propensity score overlap weight??????????????????????????????????????



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
assess_matching_balance %>% 
  ggplot(aes(x=smd, y=var, col=matched)) +
  geom_point() +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("blue","red", "green"))


################################################################################
#####                           outcome analysis                           #####
################################################################################

# separate script for this

