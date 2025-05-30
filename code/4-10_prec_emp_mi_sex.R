################################################################################

# Precarious employment and health - Understanding Society
# 4-04 - IPTW multiple imputed analytic sample - sub-group analysis by sex (FEMALE)
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
library(miceadds) # for working with imputed dataset
library(tableone) # for descriptive tables

################################################################################
#####                         load and prepare data                        #####
################################################################################


#### read in variable vectors --------------------------------------------------
source("./look_ups/variable_vectors.r")

##### load original data without imputation ------------------------------------
mi_subset2 <-  readRDS("./working_data/mi/mi_subset2.rds")

#### prepare data -------------------------------------------------------------- 

### females
df_f <- subset(mi_subset2, mi_subset2$sex_bin == 0)

df_f_str <- df_f %>% 
  summary.default() %>% as.data.frame %>% 
  dplyr::group_by(Var1) %>%  
  tidyr::spread(key = Var2, value = Freq)

# check
sapply(df_f, function(x) sum(is.na(x)))

################################################################################
#####                  Female multiple imputation model                    #####
################################################################################

#### set method type depending on variable type --------------------------------
# specify method for each incomplete variable in data
# (numeric, binary, unordered, ordered)
# ("norm", "logreg", "polyreg", "polr")
df_f_str$method <- ""
#df_f_str$method[df_f_str$Var1=="pidp"] <- ""
df_f_str$method[df_f_str$Var1=="sf12pcs_dv_t1"] <- "norm"
df_f_str$method[df_f_str$Var1=="sf12mcs_dv_t1"] <- "norm"
df_f_str$method[df_f_str$Var1=="srh_bin_t1"] <- "logreg"
df_f_str$method[df_f_str$Var1=="ghq_case4_t1"] <- "logreg"
df_f_str$method[df_f_str$Var1=="exposure1"] <- "logreg"
df_f_str$method[df_f_str$Var1=="sex_bin"] <- "pmm"
df_f_str$method[df_f_str$Var1=="age_dv_t0"] <- "norm"
df_f_str$method[df_f_str$Var1=="non_white_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="marital_status_t0"] <- "polyreg"
df_f_str$method[df_f_str$Var1=="degree_bin_t0"] <- "logreg" 
df_f_str$method[df_f_str$Var1=="gor_dv_t0"] <- "polyreg"
df_f_str$method[df_f_str$Var1=="sic2007_section_lab_t0"] <- "polyreg"
df_f_str$method[df_f_str$Var1=="soc2000_major_group_title_t0"] <- "polyreg"
df_f_str$method[df_f_str$Var1=="jbft_dv_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="small_firm_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="emp_contract_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="broken_emp_t0"] <- "polr" 
df_f_str$method[df_f_str$Var1=="j2has_dv_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="rel_pov_bin"] <- "pmm"
df_f_str$method[df_f_str$Var1=="health_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="sf12pcs_dv_t0"] <- "norm"
df_f_str$method[df_f_str$Var1=="sf12mcs_dv_t0"] <- "norm"
df_f_str$method[df_f_str$Var1=="srh_bin_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="ghq_case4_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="dep_child_bin_t0"] <- "logreg"
df_f_str$method[df_f_str$Var1=="age_dv_t1"] <- "norm" # norm if including
df_f_str$method[df_f_str$Var1=="marital_status_t1"] <- "polyreg" #polyreg if including
df_f_str$method[df_f_str$Var1=="health_t1"] <- "logreg"
df_f_str$method[df_f_str$Var1=="exp1_bin"] <- ""
df_f_str$method[df_f_str$Var1=="sex_pcs"] <- "~I(sex_bin*(sf12pcs_dv_t1-mean(sf12pcs_dv_t1)))"
df_f_str$method[df_f_str$Var1=="sex_mcs"] <- "~Isf12pcs_dv_t1-mean(sf12mcs_dv_t1, na.rm = TRUE)))"
df_f_str$method[df_f_str$Var1=="sex_srh"] <- "~I(sex_bin*srh_bin2)"
df_f_str$method[df_f_str$Var1=="sex_ghq"] <- "~I(sex_bin*ghq_bin)"
df_f_str$method[df_f_str$Var1=="age_pcs"] <- "~I(age_bin*sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
df_f_str$method[df_f_str$Var1=="age_mcs"] <- "~I(age_bin*sf12mcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
df_f_str$method[df_f_str$Var1=="age_srh"] <- "~I(age_bin*srh_bin2)"
df_f_str$method[df_f_str$Var1=="age_ghq"] <- "~I(age_bin*ghq_bin)"
df_f_str$method[df_f_str$Var1=="rel_pov_pcs"] <- "~I(rel_pov_bin*sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
df_f_str$method[df_f_str$Var1=="rel_pov_mcs"] <- "~I(rel_pov_bin*sf12mcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
df_f_str$method[df_f_str$Var1=="rel_pov_srh"] <- "~I(rel_pov_bin*srh_bin2)"
df_f_str$method[df_f_str$Var1=="rel_pov_ghq"] <- "~I(rel_pov_bin*ghq_bin)"

## check df_f_str order matches vars in df_f
sum(as.vector(df_f_str$Var1)!=as.vector(names(df_f)))

## create vector to set new default methods methods for MI
myDefaultMethod <- as.vector(df_f_str$method)

## define a custom predictorMatrix
# in matrix 1s are included in model; 0s are not
myPredictorMatrix <- make.predictorMatrix(df_f)
myPredictorMatrix[,"sex_bin"] <- 0 # add sex in as it is constant
myPredictorMatrix["sex_bin",] <- 0 # add sex in as it is constant
myPredictorMatrix[,"pidp"] <- 0
myPredictorMatrix["pidp",] <- 0
myPredictorMatrix[,"exp1_bin"] <- 0
myPredictorMatrix["exp1_bin",] <- 0
# these are needed for interaction term calculation in passive imps but should be ignored otherwise
myPredictorMatrix[,"srh_bin2"] <- 0
myPredictorMatrix["srh_bin2",] <- 0
myPredictorMatrix[,"ghq_bin"] <- 0
myPredictorMatrix["ghq_bin",] <- 0
# constituent parts of interaction terms shouldn't predict interaction term
myPredictorMatrix[c("sex_bin","sf12pcs_dv_t1"),"sex_pcs"] <- 0
myPredictorMatrix[c("sex_bin","sf12mcs_dv_t1"),"sex_pcs"] <- 0
myPredictorMatrix[c("sex_bin","srh_bin_t1"),"sex_srh"] <- 0
myPredictorMatrix[c("sex_bin","ghq_case4_t1"),"sex_ghq"] <- 0
myPredictorMatrix[c("age_dv_t0","sf12pcs_dv_t1"),"age_pcs"] <- 0
myPredictorMatrix[c("age_dv_t0","sf12mcs_dv_t1"),"age_pcs"] <- 0
myPredictorMatrix[c("age_dv_t0","srh_bin_t1"),"age_srh"] <- 0
myPredictorMatrix[c("age_dv_t0","ghq_case4_t1"),"age_ghq"] <- 0
myPredictorMatrix[c("rel_pov_bin","sf12pcs_dv_t1"),"rel_pov_pcs"] <- 0
myPredictorMatrix[c("rel_pov_bin","sf12mcs_dv_t1"),"rel_pov_pcs"] <- 0
myPredictorMatrix[c("rel_pov_bin","srh_bin_t1"),"rel_pov_srh"] <- 0
myPredictorMatrix[c("rel_pov_bin","ghq_case4_t1"),"rel_pov_ghq"] <- 0

myPredictorMatrix

#### imputation ----------------------------------------------------------------

set.seed(52267)
start_time <- Sys.time()
### for final MI model run 25 imputations and 10 iterations
## run fewer when testing code
imps2_f <- mice(
              df_f, m=25, maxit = 10, # full model
              #df_f, m=10, maxit = 5, # short model for testing
              #defaultMethod=myDefaultMethod, 
              predictorMatrix=myPredictorMatrix,
              printFlag = FALSE)
end_time <- Sys.time()
end_time - start_time

summary(imps2_f)

### check NAs are removed
## convert data from mids format to standard df
temp2 <- complete(imps2_f, "long", include = FALSE)

## sum number of NAs by var
temp2_na <- temp2 %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x)))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(df="mi")

## compare with non-imputed data
df_f_na <- df_f  %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x)))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(df="og")

temp3 <- temp2_na %>% bind_rows(df_f_na) %>% 
  pivot_wider(names_from = df, values_from = value)

write_csv(temp3, "./output/temp_output/temp3_micheck_f.csv")

## check NAs
sapply(complete(imps2_f,2), function(x) sum(is.na(x)))

#### save imputed data 
write_rds(imps2_f, "./working_data/mi/imputed_data_f.rds")



################################################################################
#####               inverse probability of treatment weighting             #####
################################################################################

imps2_f <- readRDS("./working_data/mi/imputed_data_f.rds")

#### function ------------------------------------------------------------------

iptw_func_sex <- function(data){
                   weightthem(exp1_bin ~
 #                           sex_bin +
                            age_dv_t0 +
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
                            rel_pov_bin +
                            health_t0 +
                            srh_bin_t0 +
                            ghq_case4_t0 +
                            sf12mcs_dv_t0 +
                            sf12pcs_dv_t0, 
                          data = data,
                          stabilize = TRUE,
                          estimand = "ATE",  
                          method = "ps")
}

#### females -------------------------------------------------------------------

weightit_f <- iptw_func_sex(data = imps2_f)
summary(weightit_f)

test <- bal.tab(weightit_f, un = TRUE, 
                binary = "std", continuous = "std")
test2 <- test$Balance.Across.Imputations

## probably don't need these....
#bal.plot(weightit_f, which.imp = 1, 
#         var.name = "sex_bin", 
#         which = "both")
bal.plot(weightit_f, which.imp = 1, 
         var.name = "age_dv_t0", 
         which = "both")
bal.plot(weightit_f, which.imp = 1, 
         var.name = "non_white_t0", 
         which = "both")
bal.plot(weightit_f, which.imp = 1, 
         var.name = "rel_pov_bin", 
         which = "both")

## create love plot to visualise balance between unmatched and matched data across MIs
tiff("./output/mi/sub_groups/sex/f_love_plot.tiff", width = 960)
love.plot(weightit_f, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
dev.off()
# var.names() - to clean up names
# https://www.rdocumentation.org/packages/cobalt/versions/3.4.1/topics/var.names


################################################################################
#####                               Descriptives                           #####
################################################################################

## need to re-write for this script

#### prep the data -------------------------------------------------------------

### create complete df of all imputations
weightit_f_complete <- complete(weightit_f, "long") %>% 
  # rename weights as it throws an error otherwise
  rename("ps_weights" = "weights") %>% 
  mutate(dep_child_bin_t0 = ifelse(dep_child_bin_t0==1,"Yes","No"))

### create IPTW weighted df
svy_weightit_f_complete <- svydesign(ids = ~1,
                                      data = weightit_f_complete,
                                      weights = ~ps_weights)

### create table one
table_one_f <- svyCreateTableOne(vars = cov_vector,
                               data=svy_weightit_f_complete, 
                               #  factorVars=catVars_short_vec,
                               strata = "exposure1", 
                               test =FALSE)


## non-normal numeric vector for table one
# temp df with only numeric vars
temp <- weightit_f_complete %>% dplyr::select(c(age_dv_t0,
                                                 sf12pcs_dv_t0,
                                                 sf12mcs_dv_t0))
# vector for numeric vars
nonnorm_vec <- colnames(temp)

## printed version for saving
table_one_f_sav <- print(table_one_f, showAllLevels = TRUE, smd = FALSE,
                       nonnormal = nonnorm_vec,
                       factorVars = c(catVars_vec),
                       formatOptions = list(big.mark = ",",
                                            scientific = FALSE))

### save table one
## NOTE: will have to divide number of cases by number of imputations to get pooled 
## values; other stats should be correct
write.csv(table_one_f_sav, "./output/mi/sub_groups/sex/table_one_unpooled_f.csv")

################################################################################
#####                           sort out table 1                           #####
################################################################################

table_one_f_sav <- read.csv("./output/mi/sub_groups/sex/table_one_unpooled_f.csv")

#### unexposed -----------------------------------------------------------------

## trim whitespace at both ends of strings
table_one_f_sav$unexposed <- trimws(table_one_f_sav$unexposed)

### extract first part of string
table_one_f_sav$unexposed_temp <- sub("\\ .*","",table_one_f_sav$unexposed)
## remove commas to allow conversion to numeric
table_one_f_sav$unexposed_temp <- sub(",","",table_one_f_sav$unexposed_temp)
## convert to numeric
table_one_f_sav$unexposed_temp <- as.numeric(table_one_f_sav$unexposed_temp)

## extract second part of string
table_one_f_sav$unexposed_temp2 <- sub("^[^ ]+.","",table_one_f_sav$unexposed)

## divide number of cases by number of imputations (n=25)
table_one_f_sav <- table_one_f_sav %>% 
  mutate(unexposed_temp = ifelse(X %in% c("age_dv_t0 (median [IQR])",
                                          "sf12mcs_dv_t0 (median [IQR])",
                                          "sf12pcs_dv_t0 (median [IQR])"),unexposed_temp,
                                 unexposed_temp/25))

## recreate unexposed var
table_one_f_sav$unexposed <- paste(table_one_f_sav$unexposed_temp,table_one_f_sav$unexposed_temp2)

#### exposed -------------------------------------------------------------------

## trim whitespace at both ends of strings
table_one_f_sav$exposed..employed.at.t1. <- trimws(table_one_f_sav$exposed..employed.at.t1.)

### extract first part of string
table_one_f_sav$exposed_temp <- sub("\\ .*","",table_one_f_sav$exposed..employed.at.t1.)
## remove commas to allow conversion to numeric
table_one_f_sav$exposed_temp <- sub(",","",table_one_f_sav$exposed_temp)
table_one_f_sav$exposed_temp <- sub(",","",table_one_f_sav$exposed_temp)
## convert to numeric
table_one_f_sav$exposed_temp <- as.numeric(table_one_f_sav$exposed_temp)

## extract second part of string
table_one_f_sav$exposed_temp2 <- sub("^[^ ]+.","",table_one_f_sav$exposed..employed.at.t1)

## divide number of cases by number of imputations (n=25)
table_one_f_sav <-  table_one_f_sav %>% 
  mutate(exposed_temp = ifelse(X %in% c("age_dv_t0 (median [IQR])",
                                        "sf12mcs_dv_t0 (median [IQR])",
                                        "sf12pcs_dv_t0 (median [IQR])"),exposed_temp,
                               exposed_temp/25))

## recreate exposed var
table_one_f_sav$exposed..employed.at.t1. <- paste(table_one_f_sav$exposed_temp,table_one_f_sav$exposed_temp2)

## keep only required cols
table_one_f_sav  <-  table_one_f_sav %>% dplyr::select(c("X", "level", "unexposed", "exposed..employed.at.t1."))

write.csv(table_one_f_sav, "./output/mi/sub_groups/sex/table_one_f_PAPER.csv")

################################################################################
#####                         outcome descriptives                         #####
################################################################################

#### prep the data -------------------------------------------------------------

### create complete df of all imputations
outcomes_f_complete <- weightit_f_complete %>% 
  dplyr::select(c(exposure1, 
                  sf12mcs_dv_t0, sf12mcs_dv_t1,
                  sf12pcs_dv_t0, sf12pcs_dv_t1,
                  ghq_case4_t0, ghq_case4_t1,
                  srh_bin_t0, srh_bin_t1,
                  ps_weights))

### create IPTW weighted df
svy_outcomes_f_complete <- svydesign(ids = ~1,
                                     data = outcomes_f_complete,
                                     weights = ~ps_weights)

### create table one
table_outcomes_f <- svyCreateTableOne(
                                 vars = c("sf12mcs_dv_t0", "sf12mcs_dv_t1",
                                          "sf12pcs_dv_t0", "sf12mcs_dv_t1",
                                          "ghq_case4_t0", "ghq_case4_t1",
                                          "srh_bin_t0", "srh_bin_t1"),
                                 data=svy_outcomes_f_complete, 
                                 #  factorVars=catVars_short_vec,
                                 strata = "exposure1", 
                                 test =FALSE)


## printed version for saving
table_outcomes_f_sav <- print(table_outcomes_f, showAllLevels = TRUE, smd = FALSE,
#                         nonnormal = nonnorm_vec,
#                         factorVars = c(catVars_vec),
                         formatOptions = list(big.mark = ",",
                                              scientific = FALSE))

write.csv(table_outcomes_f_sav, "./output/mi/sub_groups/sex/table_outcomes_f.csv")


################################################################################
#####                    T0 and T1 outcome descriptives                    #####
################################################################################
#
#outcomes_desc <- svyCreateTableOne(vars = c("sf12pcs_dv_t0","sf12pcs_dv_t1",
#                                            "sf12mcs_dv_t0","sf12mcs_dv_t1",
#                                            "srh_bin_t0","srh_bin_t1",
#                                            "ghq_case4_t0", "ghq_case4_t1"),
#                                   data=svy_weightit_df_complete, 
#                                   #  factorVars=catVars_short_vec,
#                                   strata = "exposure1", 
#                                   test =FALSE)
#
## printed version for saving
#outcomes_desc_sav <- print(outcomes_desc, showAllLevels = FALSE, smd = FALSE,
#                           #                      nonnormal = nonnorm_vec,
#                           factorVars = c(catVars_vec),
#                           formatOptions = list(big.mark = ",",
#                                                scientific = FALSE))
#
#write.csv(outcomes_desc_sav, "./output/mi/weighted_descriptives/outcomes_desc.csv")
#

################################################################################
#####                      IPTW double-robust model                        #####
################################################################################

#### functions -----------------------------------------------------------------
### sf-12 pcs
iptw_dr_pcs <- function(data){
                      with(data,
                          glmmTMB(sf12pcs_dv_t1 ~
                                    exposure1 +
                                    sf12pcs_dv_t0 +
#                                    sex_bin +
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
                                    rel_pov_bin +
                                    health_t0 +
                                    health_t1 +
                                    # interaction terms
#                                    sex_bin*age_dv_t0 +
#                                    sex_bin*rel_pov_bin +
                                    age_dv_t0*rel_pov_bin +
                                    (1|pidp)))
}

### sf-12 mcs
iptw_dr_mcs <- function(data){
  with(data,
       glmmTMB(sf12mcs_dv_t1 ~
                 exposure1 +
                 sf12mcs_dv_t0 +
                 #                                    sex_bin +
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
                 rel_pov_bin +
                 health_t0 +
                 health_t1 +
                 # interaction terms
                 #                                    sex_bin*age_dv_t0 +
                 #                                    sex_bin*rel_pov_bin +
                 age_dv_t0*rel_pov_bin +
                 (1|pidp)))
}

### self-rated health
iptw_dr_srh <- function(data){
  with(data,
       glmmTMB(srh_bin_t1 ~
                 exposure1 +
                 srh_bin_t0 +
                 #                                    sex_bin +
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
                 rel_pov_bin +
                 health_t0 +
                 health_t1 +
                 # interaction terms
                 #                                    sex_bin*age_dv_t0 +
                 #                                    sex_bin*rel_pov_bin +
                 age_dv_t0*rel_pov_bin +
                 (1|pidp),
               family=binomial(link="logit")))
}

### ghq-12
iptw_dr_ghq <- function(data){
  with(data,
       glmmTMB(ghq_case4_t1 ~
                 exposure1 +
                 ghq_case4_t0 +
                 #                                    sex_bin +
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
                 rel_pov_bin +
                 health_t0 +
                 health_t1 +
                 # interaction terms
                 #                                    sex_bin*age_dv_t0 +
                 #                                    sex_bin*rel_pov_bin +
                 age_dv_t0*rel_pov_bin +
                 (1|pidp),
                 family=binomial(link="logit")))
}

#### females -------------------------------------------------------------------

### SF-12 PCS --------------
f_mods_pcs <- iptw_dr_pcs(data=weightit_f)


f_pooled_pcs <- pool(f_mods_pcs)

weightit_pooled_pcs_df <- data.frame(summary(f_pooled_pcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)


#### SF-12 MCS -----------------------------------------------------------------

f_mods_mcs <- iptw_dr_mcs(data=weightit_f)


f_pooled_mcs <- pool(f_mods_mcs)

weightit_pooled_mcs_df <- data.frame(summary(f_pooled_mcs, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)


#### Self-rated health ---------------------------------------------------------

f_mods_srh <- iptw_dr_srh(data=weightit_f)


f_pooled_srh <- pool(f_mods_srh)

weightit_pooled_srh_df <- data.frame(summary(f_pooled_srh, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)


#### GHQ-12 caseness -----------------------------------------------------------

f_mods_ghq <- iptw_dr_ghq(data=weightit_f)


f_pooled_ghq <- pool(f_mods_ghq)

weightit_pooled_ghq_df <- data.frame(summary(f_pooled_ghq, conf.int = TRUE)) %>% 
  rename(lci = X2.5..,
         uci = X97.5..)


################################################################################
#####                   combine estimates into single df                   #####
################################################################################

pcs_df <- weightit_pooled_pcs_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "SF-12 PCS",
         est_type = "coefficient")

mcs_df <- weightit_pooled_mcs_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "SF-12 MCS",
         est_type = "coefficient")

srh_df <- weightit_pooled_srh_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "Poor self-rated health",
         est_type = "OR")

ghq_df <- weightit_pooled_ghq_df %>% 
  filter(term == "exposure1exposed (employed at t1)") %>% 
  dplyr::select(-c(component, statistic,df)) %>% 
  mutate(outcome = "GHQ-12 caseness (4+)",
         est_type = "OR")

combined_df <- pcs_df %>% 
  bind_rows(mcs_df, srh_df, ghq_df) %>% 
  dplyr::select(c(outcome,	est_type,	estimate,	std.error,	p.value,	lci,	uci)) %>% 
  mutate(sub_group = "female")

write.csv(combined_df, "./output/mi/weighted_outcomes/mi_dr_iptw_df_sex_f.csv")


##################


sapply(complete(imps2_f,"long"), function(x) sum(is.na(x)))
