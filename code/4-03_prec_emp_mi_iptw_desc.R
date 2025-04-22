################################################################################

# Precarious employment and health - Understanding Society
# 4-03 -  Descriptive analysis for IPTW multiple imputed outcome models 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a)  produces love plot to show balance between treatment groups
# (b)  produces Table 1 descriptives


################################################################################

## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
library(tableone) # for creating table 1
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

##### load original data without imputation ------------------------------------
#mi_subset2 <-  readRDS("./working_data/mi/mi_subset2.rds")

#### load imputed data --------------------------------------------------------
weightit_df <- readRDS("./working_data/mi/weightit_df.rds")

## check forNAs in imputed data
sapply(complete(weightit_df,"long"), function(x) sum(is.na(x)))

sapply(complete(weightit_df,"long"), function(x) sum(x=="missing"))


missing_check <- complete(weightit_df,"long")

#### prepare data -------------------------------------------------------------- 


################################################################################
#####                 check balance between treatment groups               #####
################################################################################


test <- bal.tab(weightit_df, un = TRUE, 
                binary = "std", continuous = "std")
test2 <- test$Balance.Across.Imputations

## probably don't need these....
bal.plot(weightit_df, which.imp = 1, 
         var.name = "sex_bin", 
         which = "both")
bal.plot(weightit_df, which.imp = 1, 
         var.name = "age_dv_t0", 
         which = "both")
bal.plot(weightit_df, which.imp = 1, 
         var.name = "non_white_t0", 
         which = "both")

## create love plot to visualise balance between unmatched and matched data across MIs
tiff("./output/mi/mi_descriptives/iptw_love_plot.tiff", width = 980)
love.plot(weightit_df, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
# var.names() - to clean up names
dev.off()

pdf("./output/mi/mi_descriptives/iptw_love_plot.pdf", width = 980, height = 980)
love.plot(weightit_df, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
# var.names() - to clean up names
dev.off()

################################################################################
#####                               Descriptives                           #####
################################################################################

#### prep the data -------------------------------------------------------------

### create complete df of all imputations
weightit_df_complete <- complete(weightit_df, "long") %>% 
# rename weights as it throws an error otherwise
    rename("ps_weights" = "weights") %>% 
  mutate(dep_child_bin_t0 = ifelse(dep_child_bin_t0==1,"Yes","No"))

### create IPTW weighted df
svy_weightit_df_complete <- svydesign(ids = ~1,
                      data = weightit_df_complete,
                      weights = ~ps_weights)

### create table one
table_one <- svyCreateTableOne(vars = c(cov_vector),
  data=svy_weightit_df_complete, 
#  factorVars=catVars_short_vec,
  strata = "exposure1", 
  test =FALSE)


## non-normal numeric vector for table one
# temp df with only numeric vars
temp <- weightit_df_complete %>% dplyr::select(c(age_dv_t0,
                                             sf12pcs_dv_t0,
                                             sf12mcs_dv_t0))
# vector for numeric vars
nonnorm_vec <- colnames(temp)

## printed version for saving
table_one_sav <- print(table_one, showAllLevels = TRUE, smd = FALSE,
                       nonnormal = "age_dv_t0",
                       factorVars = c(catVars_vec),
                       formatOptions = list(big.mark = ",",
                                            scientific = FALSE))

### save table one
## NOTE: will have to divide number of cases by number of imputations to get pooled 
## values; other stats should be correct
write.csv(table_one_sav, "./output/mi/weighted_descriptives/table_one_unpooled.csv")


################################################################################
#####                    T0 and T1 outcome descriptives                    #####
################################################################################

#### not using
nonnorm_vec2 <- c("sf12pcs_dv_t0","sf12pcs_dv_t1",
                                            "sf12mcs_dv_t0","sf12mcs_dv_t1")

## create table
outcomes_desc <- svyCreateTableOne(vars = c("sf12pcs_dv_t0","sf12pcs_dv_t1",
                                            "sf12mcs_dv_t0","sf12mcs_dv_t1",
                                            "srh_bin_t0","srh_bin_t1",
                                            "ghq_case4_t0", "ghq_case4_t1"),
                                               data=svy_weightit_df_complete, 
                                   #  factorVars=catVars_short_vec,
                                               strata = "exposure1", 
                                               test =FALSE)

## printed version for saving
outcomes_desc_sav <- print(outcomes_desc, showAllLevels = FALSE, smd = FALSE,
#                           nonnormal = nonnorm_vec2,
#                           factorVars = c(catVars_vec),
                       formatOptions = list(big.mark = ",",
                                            scientific = FALSE))

write.csv(outcomes_desc_sav, "./output/mi/weighted_descriptives/outcomes_desc.csv")

################################################################################
#####                           sort out table 1                           #####
################################################################################

table_one_sav <- read.csv("./output/mi/weighted_descriptives/table_one_unpooled.csv")

#### unexposed -----------------------------------------------------------------

## trim whitespace at both ends of strings
table_one_sav$unexposed <- trimws(table_one_sav$unexposed)

### extract first part of string
table_one_sav$unexposed_temp <- sub("\\ .*","",table_one_sav$unexposed)
## remove commas to allow conversion to numeric
table_one_sav$unexposed_temp <- sub(",","",table_one_sav$unexposed_temp)
## convert to numeric
table_one_sav$unexposed_temp <- as.numeric(table_one_sav$unexposed_temp)

## extract second part of string
table_one_sav$unexposed_temp2 <- sub("^[^ ]+.","",table_one_sav$unexposed)

## divide number of cases by number of imputations (n=25)
table_one_sav <- table_one_sav %>% 
  mutate(unexposed_temp = ifelse(X %in% c("age_dv_t0 (median [IQR])",
                                          "sf12mcs_dv_t0 (mean (SD))",
                                          "sf12pcs_dv_t0 (mean (SD))"),unexposed_temp,
                                 unexposed_temp/25))

## recreate unexposed var
table_one_sav$unexposed <- paste(table_one_sav$unexposed_temp,table_one_sav$unexposed_temp2)

#### exposed -------------------------------------------------------------------

## trim whitespace at both ends of strings
table_one_sav$exposed..employed.at.t1. <- trimws(table_one_sav$exposed..employed.at.t1.)

### extract first part of string
table_one_sav$exposed_temp <- sub("\\ .*","",table_one_sav$exposed..employed.at.t1.)
## remove commas to allow conversion to numeric
table_one_sav$exposed_temp <- sub(",","",table_one_sav$exposed_temp)
table_one_sav$exposed_temp <- sub(",","",table_one_sav$exposed_temp)
## convert to numeric
table_one_sav$exposed_temp <- as.numeric(table_one_sav$exposed_temp)

## extract second part of string
table_one_sav$exposed_temp2 <- sub("^[^ ]+.","",table_one_sav$exposed..employed.at.t1)

## divide number of cases by number of imputations (n=25)
table_one_sav <-  table_one_sav %>% 
  mutate(exposed_temp = ifelse(X %in% c("age_dv_t0 (median [IQR])",
                                          "sf12mcs_dv_t0 (mean (SD))",
                                          "sf12pcs_dv_t0 (mean (SD))"),exposed_temp,
                                 exposed_temp/25))

## recreate exposed var
table_one_sav$exposed..employed.at.t1. <- paste(table_one_sav$exposed_temp,table_one_sav$exposed_temp2)

## keep only required cols
table_one_sav  <-  table_one_sav %>% dplyr::select(c("X", "level", "unexposed", "exposed..employed.at.t1."))

write.csv(table_one_sav, "./output/mi/weighted_descriptives/table_one_PAPER.csv")
