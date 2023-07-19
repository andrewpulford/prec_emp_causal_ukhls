################################################################################

# Precarious employment and health - Understanding Society
# 3-02 - Paired unweighted descriptives 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) produces unweighted descriptive analysis for complete case paired data


################################################################################

## remove any existing objects from global environment
rm(list=ls()) 

## turn off scientific notation
options(scipen = 999)

################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
library(Hmisc) # histogram plotting across cols
library(tableone) # for creating table one

################################################################################
#####                         load and prepare data                        #####
################################################################################

####load analytic df -----------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")

# check no NAs
sapply(pair_cc_analytic, function(x) sum(is.na(x)))

#pair_cc_analytic$fimnnet_dv_t0 <- as.numeric(as.character(pair_cc_analytic$fimnnet_dv_t0))


#pair_cc_analytic$srh_dv_t0 <- factor(pair_cc_analytic$srh_dv_t0, 
#                                     levels = c("excellent", "very good", "good",
#                                                "fair", "poor"))

pair_cc_analytic$sf12pcs_dv_t0 <- as.character(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t0)

pair_cc_analytic$sf12mcs_dv_t0 <- as.character(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t0)


#####----------------------------------------------------------------------#####
#####                          Save sample chars data                      #####
#####----------------------------------------------------------------------#####

## as dataframe
#write_rds(sample_chars, "./working_data/sample_chars.rds")
#
#
#
##### job loss between t0-t1
#
#emp_status_t1 <- pair_cc_analytic %>%
#  group_by(exposure1)  %>%  
#  summarise(n=n()) %>%  
#  mutate(est = n/sum(n)*100) %>% 
#  mutate(var="Employment status t1") %>% 
#  rename("measure"= "exposure1") %>% 
#  dplyr::select(var, measure, n, est)
#
## if unemployed at t0 or if unemp spell between t0-1
#lost_job <- pair_cc_analytic %>% 
#  group_by(exposure2) %>% 
#  summarise(n=n()) %>%  
#  mutate(est = n/sum(n)*100) %>% 
#  mutate(var="Job loss between t0 and t1") %>% 
#  rename("measure"= "exposure2") %>% 
#  dplyr::select(var, measure, n, est)
#
#
################################################################################
#####        Table 1: participant characteristics by exposure group        #####
################################################################################

pair_cc_analytic <- droplevels(pair_cc_analytic)

#pair_cc_analytic <- pair_cc_analytic %>% 
#  mutate(across(.cols = everything(), 
#                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 

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

pair_cc_analytic <- pair_cc_analytic %>% 
  dplyr::select(all_of(cov_vector), exposure1, exposure2)

## histogram for each numeric variable
# temp df with only numeric vars
temp <- pair_cc_analytic %>% dplyr::select(where(is.numeric))

# histogram
#hist.data.frame(temp)

# normality test for age
temp2 <- sample (temp$age_dv_t0, size=5000)
shapiro.test(temp2)

# winsorized version of income to deal with extreme values
#hist(DescTools::Winsorize(pair_cc_analytic$fimnnet_dv_t0))

## vector for numeric vars
nonnorm_vec <- colnames(temp)

## unemployed at T1
table_one <- tableone::CreateTableOne(vars = cov_vector, strata = "exposure1",
                            data = pair_cc_analytic,
                            test = FALSE)

table_one_sav <- print(table_one, showAllLevels = TRUE, smd = TRUE, 
                       nonnormal = nonnorm_vec,
                       formatOptions = list(big.mark = ","))

# Count covariates with important imbalance
table_one_smd <- data.frame(ExtractSmd(table_one))
table_one_smd <- table_one_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_smd, "./working_data/table_one_unmatched_smd.csv")

addmargins(table(ExtractSmd(table_one) > 0.1))

# plot unmatched smd's
table_one_smd %>% 
  ggplot(aes(x=smd, y=var, col=imbalance_flag)) +
  geom_point() +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("blue","red"))

## job loss between t0 and t1
table_one_alt <- CreateTableOne(vars = cov_vector, strata = "exposure2", 
                                data = pair_cc_analytic)

# Count covariates with important imbalance
table_one_alt_smd <- data.frame(ExtractSmd(table_one_alt))
table_one_alt_smd <- table_one_alt_smd %>% 
  rownames_to_column("var") %>% # Apply rownames_to_column
  rename("smd" = "X1.vs.2") %>% 
  mutate(imbalance_flag = ifelse(smd>0.1,"SMD>0.1","SMD<=0.1"),
         matched = "unmatched")

write.csv(table_one_alt_smd, "./output/unmatched_descriptives/table_one_alt_unmatched_smd.csv")

# plot unmatched smd's
table_one_alt_smd %>% 
  ggplot(aes(x=smd, y=var, col=imbalance_flag)) +
  geom_point() +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  theme_bw() +
  scale_color_manual(values = c("blue","red"))

addmargins(table(ExtractSmd(table_one_alt) > 0.1))

table_one_alt_sav <- print(table_one_alt, showAllLevels = TRUE, smd = TRUE, 
                           formatOptions = list(big.mark = ","))

### save tables

write.csv(table_one_sav, "./output/unmatched_descriptives/table_one_unmatched.csv")
write.csv(table_one_alt_sav, "./output/unmatched_descriptives/table_one_alt_unmatched.csv")

