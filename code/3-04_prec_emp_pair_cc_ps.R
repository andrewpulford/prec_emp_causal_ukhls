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
library(glmmTMB) # for multi-level modelling (faster than lme4)
library(numDeriv) # for gradient checks

################################################################################
#####                         load and prepare data                        #####
################################################################################

####load eligible cases --------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")

### convert binary outcome vars to factors to allow svyglm to work
pair_cc_analytic$srh_bin_t1 <- factor(pair_cc_analytic$srh_bin_t1)
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
                "rel_pov_t0",
                "health_t0",
                "srh_bin_t0",
                "ghq_case4_t0",
                "sf12mcs_dv_t0",
                "sf12pcs_dv_t0")

cov_vector2 <- c("age_dv_t1",  
                "marital_status_t1",
                "gor_dv_t1",
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

nonnorm_vec <- (c("age_dv_t0", "sf12mcs_dv_t0", "sf12pcs_dv_t0"))

#### keep only variables required for propensity score (and other analysis) ----
pair_cc_ps <- pair_cc_analytic %>% 
  dplyr::select(pidp, all_of(c(cov_vector, cov_vector2, outcome_vector))) #%>% 
#  dplyr::select(-c(psu, strata, wt_name, wt_value))

#### convert pidp to factor 
pair_cc_ps$pidp <- factor(pair_cc_ps$pidp)

#### convert scale vars to numeric ---------------------------------------------
pair_cc_ps <- pair_cc_ps %>% 
#  mutate(fimnnet_dv_t0 = as.numeric(fimnnet_dv_t0)) %>%
  ## these to character first so score converted rather than factor level
  mutate(sf12mcs_dv_t0 = as.character(sf12mcs_dv_t0),
         sf12pcs_dv_t0 = as.character(sf12pcs_dv_t0),
         sf12mcs_dv_t1 = as.character(sf12mcs_dv_t1),
         sf12pcs_dv_t1 = as.character(sf12pcs_dv_t1)) %>% 
  mutate(sf12mcs_dv_t0 = as.numeric(sf12mcs_dv_t0),
         sf12pcs_dv_t0 = as.numeric(sf12pcs_dv_t0),
         sf12mcs_dv_t1 = as.numeric(sf12mcs_dv_t1),
         sf12pcs_dv_t1 = as.numeric(sf12pcs_dv_t1))

## change any character vars to factor
pair_cc_ps[sapply(pair_cc_ps, is.character)] <- lapply(pair_cc_ps[sapply(pair_cc_ps, is.character)], as.factor)

#### revel exposure vars so that unexposed is lower ref category ---------------
pair_cc_ps$exposure1 <- factor(pair_cc_ps$exposure1,
                                     levels = rev(levels(pair_cc_ps$exposure1)))

pair_cc_ps$exposure2 <- factor(pair_cc_ps$exposure2,
                                     levels = rev(levels(pair_cc_ps$exposure2)))


#### remove any NAs not coded as a missing category ----------------------------
##  check NAs
#sapply(pair_cc_ps, function(x) sum(is.na(x)))
## remove
#pair_cc_ps <- na.omit(pair_cc_ps)
## check again
#sapply(pair_cc_ps, function(x) sum(is.na(x)))


### use this if full df is too big
#test_df <- pair_cc_ps %>% 
#  slice_sample(prop = 1)

#### load unmatched SMD df
table_one_unmatched_smd <- read.csv("./working_data/table_one_unmatched_smd.csv") %>% 
  dplyr::select(c(var, smd, imbalance_flag, matched))


################################################################################
#####                               functions                              #####
################################################################################

#### propensity score model - basic --------------------------------------------
ps_model <- function(data = pair_cc_ps, outcome){
  glm(outcome ~
          sex_dv_t0 +
          age_dv_t0 +
          non_white_t0  +
          marital_status_t0 +
          hiqual_dv_t0 +
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
          srh_bin_t0 +
          ghq_case4_t0 +
          sf12mcs_dv_t0 +
          sf12pcs_dv_t0,
        family = binomial(link="logit"),
        data = data)
  
  
}

#### propensity score model - mlm ----------------------------------------------
ps_model_mlm <- function(data = pair_cc_ps, outcome){
  glmmTMB::glmmTMB(outcome ~
         sex_dv_t0 +
         age_dv_t0 +
         non_white_t0  +                
         marital_status_t0 +            
         hiqual_dv_t0 +                 
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
         srh_bin_t0 +                   
         ghq_case4_t0 +                 
         sf12mcs_dv_t0 +                
         sf12pcs_dv_t0 +                
         (1|pidp),
            family = binomial(link="logit"),
            data = data)#, control=glmerControl(optimizer="bobyqa",
#                     optCtrl=list(maxfun=2e5)))
  
  
}




################################################################################
#####                     Unemployment at t1 PS model                      #####
################################################################################

### call the function (and benchmark time - takes a while to run foe MLM ~2 days!)
## not MLM
start_time <- Sys.time()
ps_mod_exp1_noMLM <- ps_model(data = pair_cc_ps, outcome = pair_cc_ps$exposure1)
end_time <- Sys.time()
end_time - start_time

summary(ps_mod_exp1_noMLM)

## MLM
start_time <- Sys.time()
ps_mod_exp1 <- ps_model_mlm(data = pair_cc_ps, outcome = pair_cc_ps$exposure1)
end_time <- Sys.time()
end_time - start_time

### summary of model
summary(ps_mod_exp1)

### checks
#
## rescale vars
#numcols <- grep("^c\\.",names(pair_cc_ps))
#pair_cc_ps_scaled <- pair_cc_ps
#pair_cc_ps_scaled[,numcols] <- scale(pair_cc_ps_scaled[,numcols])
#ps_mod_exp1 <- update(ps_mod_exp1,data=pair_cc_ps_scaled)#
#
## Check singularity
#tt <- getME(ps_mod_exp1,"theta")
#ll <- getME(ps_mod_exp1,"lower")
#min(tt[ll==0])
#
## Double-checking gradient calculations
#derivs1 <- ps_mod_exp1@optinfo$derivs
#check_grad1 <- with(derivs1,solve(Hessian,gradient))
#max(abs(check_grad1))
#
# try with numDeriv
#dd <- update(ps_mod_exp1,devFunOnly=TRUE)
#pars <- unlist(getME(ps_mod_exp1,c("theta","fixef")))
#grad2 <- grad(dd,pars)
#hess2 <- hessian(dd,pars)
#check_grad2 <- solve(hess2,grad2)
#max(pmin(abs(check_grad2),abs(grad2)))

## try update with max iterations
#ss <- getME(ps_mod_exp1,c("theta","fixef"))
#start_time <- Sys.time()
#ps_mod_exp1_update <- update(ps_mod_exp1,start=ss,
#                             control=glmerControl(optimizer="bobyqa",
#                                                  optCtrl=list(maxfun=2e5)))
#end_time <- Sys.time()
#end_time - start_time



## check for 1 or 0 predicted probabilities
pair_cc_ps$y_pred <- predict(ps_mod_exp1, pair_cc_ps, type="response")

summary(pair_cc_ps$y_pred)

test <- pair_cc_ps %>% filter(y_pred==min(y_pred)) %>% 
  dplyr::select(pidp, y_pred)

### predicted probability of being assigned to exposed group
pair_cc_ps$ps_exp1 <- predict(ps_mod_exp1, type = "response")
summary(pair_cc_ps$ps_exp1)

### predicted probability of being assigned to unexposed group
pair_cc_ps$ps_noexp1 <- 1-pair_cc_ps$ps_exp1
summary(pair_cc_ps$ps_noexp1)

### predicted probability of being assigned to actual exposure status
pair_cc_ps$ps_assign <- NA
pair_cc_ps$ps_assign[pair_cc_ps$exposure1=="exposed (employed at t1)"] <- pair_cc_ps$ps_exp1[pair_cc_ps$exposure1=="exposed (employed at t1)"]
pair_cc_ps$ps_assign[pair_cc_ps$exposure1=="unexposed"] <- pair_cc_ps$ps_noexp1[pair_cc_ps$exposure1=="unexposed"]

### smaller of ps_exp1 and ps_noexp1 for matchnig weight

pair_cc_ps$ps_min <- pmin(pair_cc_ps$ps_exp1, pair_cc_ps$ps_noexp1)

################################################################################
#####                       propensity score matching                      #####
################################################################################

### remove missing as category for tables 
#pair_cc_ps <- pair_cc_ps %>% 
#  mutate(across(.cols = everything(), 
#                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 


#### match propensity scores using Matching package ----------------------------
list_match <- Match(Tr = (pair_cc_ps$exposure1=="exposed (employed at t1)"),
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

table_one_weighted_sav <- print(table_one_weighted, showAllLevels = TRUE, smd = TRUE,
                                nonnormal = nonnorm_vec,
                                formatOptions = list(big.mark = ","))

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

################################################################################
#IPW package test 
#
#library(ipw)
#
#### remove missing as category for tables 
#pair_cc_ps <- pair_cc_ps %>% 
#  mutate(across(.cols = everything(), 
#                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 
#
#pair_cc_ipw <- pair_cc_ps %>% 
#  mutate(exposure1 = as.numeric(exposure1),
#         sex_dv_t0 = as.numeric(sex_dv_t0)) %>% 
#  mutate(exposure1 = exposure1-1,
#         male_t0 = sex_dv_t0-1,
#         non_white_t0 = ifelse(non_white_t0==3,0,1)) %>% 
#  drop_na(everything())
##  dplyr::select(male_t0) %>% 
##  group_by(male_t0) %>% 
##  summarise(n=n())
#
#weight <- ipwpoint(exposure = exposure1, family = "binomial", link = "logit",
#                   numerator =~ 1,
#                   denominator =~ male_t0 + age_dv_t0 + non_white_t0,
#                   id = pair_cc_ipw$pidp,
#                   trunc = 0.01, data = as.data.frame(pair_cc_ipw))
##currentDataset$.ipw0 = weight$weights.trunc
#
#
#
#  marital_status_t0 +
#  hiqual_dv_t0 +
#  gor_dv_t0 +
#  sic2007_section_lab_t0 +
#  soc2000_major_group_title_t0 +
#  jbft_dv_t0 +
#  small_firm_t0 +
#  #         jbhrs_t0 +
#  emp_contract_t0 +
#  broken_emp_t0 +
#  j2has_dv_t0 +
#  rel_pov_t0 +
#  health_t0 +
#  srh_bin_t0 +
#  ghq_case4_t0 +
#  sf12mcs_dv_t0 +
#  sf12pcs_dv_t0#