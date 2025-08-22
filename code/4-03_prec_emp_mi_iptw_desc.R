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
love.plot(weightit_df, thresholds = 0.1, 
          var.names = c(sex_bin = "Gender",
            age_dv_t0 = "Age (Years)",
            `non_white_t0_Non-white` = "Non-white ethnicity",
            `marital_status_t0_single` = "Marital status: Single",
            `marital_status_t0_divorced/separated/widowed` = "Marital status: Divorced/separated/widowed",
            dep_child_bin_t0 = "Has dependent children",
            `degree_bin_t0_Degree or higher` = "Highest qualification: Degree or higher",
            `gor_dv_t0_east midlands` = "Region: East Midlands",
            `gor_dv_t0_east of england` = "Region: East of England",
            `gor_dv_t0_london` = "Region: London",
            `gor_dv_t0_north east` = "Region: North East",
            `gor_dv_t0_north west` = "Region: North West",
            `gor_dv_t0_northern ireland` = "Region: Northern Ireland",
            `gor_dv_t0_scotland` = "Region: Scotland",
            `gor_dv_t0_south east` = "Region: South East",
            `gor_dv_t0_south west` = "Region: South West",
            `gor_dv_t0_wales` = "Region: Wales",
            `gor_dv_t0_west midlands` = "Region: West Midlands",
            `gor_dv_t0_yorkshire and the humber` = "Region: Yorkshire and the Humber",
            `sic2007_section_lab_t0_Accommodation and food service activities` = "Industry: Accommodation and food service activities",
            `sic2007_section_lab_t0_Administrative and support service activities` = "Industry: Administrative and support service activities",
            `sic2007_section_lab_t0_Construction` = "Industry: Construction",
            `sic2007_section_lab_t0_Education` = "Industry: Education",
            `sic2007_section_lab_t0_Human health and social work activities` = "Industry: Human health and social work activities",
            `sic2007_section_lab_t0_Manufacturing` = "Industry: Manufacturing",
            `sic2007_section_lab_t0_Other industry` = "Industry: Other industry",
            `sic2007_section_lab_t0_Professional, scientific and technical activities` = "Industry: Professional, scientific and technical activities",
            `sic2007_section_lab_t0_Public administration and defence; compulsory social security` = "Industry: Public administration and defence; compulsory social security",
            `sic2007_section_lab_t0_Transportation and storage` = "Industry: Transportation and storage",
            `sic2007_section_lab_t0_` = "Industry: Wholesale and retail trade; repair of motor vehicles and motorcycles",
            `soc2000_major_group_title_t0_Administrative and secretarial occupations` = "Occupation: Administrative and secretarial",
            `soc2000_major_group_title_t0_Associate professional and technical occupations` = "Occupation: Associate professional and technical",
            `soc2000_major_group_title_t0_Elementary occupations` = "Occupation: Elementary",
            `soc2000_major_group_title_t0_Managers and senior officials` = "Occupation: Managers and senior officials",
            `soc2000_major_group_title_t0_Personal service occupations` = "Occupation: Personal service",
            `soc2000_major_group_title_t0_Process, plant and machine operatives` = "Occupation: Process, plant and machine operatives",
            `soc2000_major_group_title_t0_Sales and customer service occupations` = "Occupation: Sales and customer service",
            `soc2000_major_group_title_t0_Science and technology professionals` = "Occupation: Science and technology professionals",
            `soc2000_major_group_title_t0_Skilled trade occupations` = "Occupation: Skilled trade",
            `jbft_dv_t0_PT employee` = "Part-time worker",
            `small_firm_t0_under 50 employees` = "Works for small firm (<50 employees)",
            `emp_contract_t0_fixed-term` = "Fixed-term employment contract",
            `broken_emp_t0_Broken employment` = "Emplyment discontinuity",
            j2has_t0_yes = "In multiple employment",
            rel_pov_bin = "In relative poverty",
            health_t0_yes = "Long-standing illness or impairment",
            `srh_bin_t0_fair/poor` = "Fair/poor self-rated health",
            `ghq_case4_t0_4 or more` = "GHQ-12 caseness",
            sf12mcs_dv_t0 = "SF-12 Mental Component Summary (MCS)",
            sf12pcs_dv_t0 = "SF-12 Physical Component Summary (PCS)"),
          stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
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
  mutate(dep_child_bin_t0 = ifelse(dep_child_bin_t0==1,"Yes","No"),
         sex_bin = ifelse(sex_bin==0, "Female","Male"),
         rel_pov_bin = ifelse(rel_pov_bin==0,"No","Yes"))

### create IPTW weighted df
svy_weightit_df_complete <- svydesign(ids = ~1,
                      data = weightit_df_complete,
                      weights = ~ps_weights)

### create table one
table_one <- svyCreateTableOne(vars = c("sex_bin", cov_vector, "rel_pov_bin"),
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
                       factorVars = c(catVars_vec, "sex_bin", "rel_pov_bin"),
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
                                            "srh_bin_t0","srh_bin2",
                                            "ghq_case4_t0", "ghq_bin"),
                                               data=svy_weightit_df_complete, 
                                   factorVars = c("srh_bin_t0", "srh_bin2", "ghq_bin"),
                                   strata = "exposure1", 
                                               test =FALSE)

## printed version for saving
outcomes_desc_sav <- print(outcomes_desc, showAllLevels = FALSE, smd = FALSE,
#                           nonnormal = nonnorm_vec2,
                           factorVars = c("srh_bin_t0", "srh_bin2", "ghq_bin"),
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

################################################################################
#####                        sort out outcomes table                       #####
################################################################################

outcomes_desc_sav <- read.csv("./output/mi/weighted_descriptives/outcomes_desc.csv")

#### unexposed -----------------------------------------------------------------

## trim whitespace at both ends of strings
outcomes_desc_sav$unexposed <- trimws(outcomes_desc_sav$unexposed)

### extract first part of string
outcomes_desc_sav$unexposed_temp <- sub("\\ .*","",outcomes_desc_sav$unexposed)
## remove commas to allow conversion to numeric
outcomes_desc_sav$unexposed_temp <- sub(",","",outcomes_desc_sav$unexposed_temp)
## convert to numeric
outcomes_desc_sav$unexposed_temp <- as.numeric(outcomes_desc_sav$unexposed_temp)

## extract second part of string
outcomes_desc_sav$unexposed_temp2 <- sub("^[^ ]+.","",outcomes_desc_sav$unexposed)

## divide number of cases by number of imputations (n=25)
outcomes_desc_sav <- outcomes_desc_sav %>% 
  mutate(unexposed_temp = ifelse(X %in% c("sf12mcs_dv_t0 (mean (SD))",
                                          "sf12mcs_dv_t1 (mean (SD))",
                                          "sf12pcs_dv_t0 (mean (SD))",
                                          "sf12pcs_dv_t1 (mean (SD))"),unexposed_temp,
                                 unexposed_temp/25))

## recreate unexposed var
outcomes_desc_sav$unexposed <- paste(outcomes_desc_sav$unexposed_temp,outcomes_desc_sav$unexposed_temp2)

#### exposed -------------------------------------------------------------------

## trim whitespace at both ends of strings
outcomes_desc_sav$exposed..employed.at.t1. <- trimws(outcomes_desc_sav$exposed..employed.at.t1.)

### extract first part of string
outcomes_desc_sav$exposed_temp <- sub("\\ .*","",outcomes_desc_sav$exposed..employed.at.t1.)
## remove commas to allow conversion to numeric
outcomes_desc_sav$exposed_temp <- sub(",","",outcomes_desc_sav$exposed_temp)
outcomes_desc_sav$exposed_temp <- sub(",","",outcomes_desc_sav$exposed_temp)
## convert to numeric
outcomes_desc_sav$exposed_temp <- as.numeric(outcomes_desc_sav$exposed_temp)

## extract second part of string
outcomes_desc_sav$exposed_temp2 <- sub("^[^ ]+.","",outcomes_desc_sav$exposed..employed.at.t1)

## divide number of cases by number of imputations (n=25)
outcomes_desc_sav <-  outcomes_desc_sav %>% 
  mutate(exposed_temp = ifelse(X %in% c("sf12mcs_dv_t0 (mean (SD))",
                                          "sf12mcs_dv_t1 (mean (SD))",
                                          "sf12pcs_dv_t0 (mean (SD))",
                                          "sf12pcs_dv_t1 (mean (SD))"),exposed_temp,
                                 exposed_temp/25))

## recreate exposed var
outcomes_desc_sav$exposed..employed.at.t1. <- paste(outcomes_desc_sav$exposed_temp,outcomes_desc_sav$exposed_temp2)

## keep only required cols
outcomes_desc_sav  <-  outcomes_desc_sav %>% dplyr::select(c("X", "unexposed", "exposed..employed.at.t1."))

write.csv(outcomes_desc_sav, "./output/mi/weighted_descriptives/outcomes_desc_PAPER.csv")
