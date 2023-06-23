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
library(tableone) # for creating table one

################################################################################
#####                         load and prepare data                        #####
################################################################################

####load analytic df -----------------------------------------------------------
pair_cc_analytic <- readRDS("./working_data/pair_cc_analytic.rds")


################################################################################
##### descriptives ######
################################################################################

#####----------------------------------------------------------------------#####
#####                     Personal characteristics                         #####
#####----------------------------------------------------------------------#####

#### sex -----------------------------------------------------------------------
#pair_cc_analytic$sex_dv <- droplevels(pair_cc_analytic$sex_dv)

sex <- pair_cc_analytic %>% group_by(sex_dv_t0) %>% summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Sex") %>% 
  rename("measure"= "sex_dv_t0") %>% 
  dplyr::select(var, measure, n, est) %>% 
  arrange(factor(measure, levels = c("Female","Male")))

sample_chars <- sex

#### age -----------------------------------------------------------------------
age_mean <- pair_cc_analytic %>%  
  summarise(est = mean(as.numeric(as.numeric(age_dv_t0)), na.rm = TRUE)) %>% 
  mutate(var="Age", measure="Mean", n=NA) %>% 
  dplyr::select(var, measure, n, est)

sample_chars <- sample_chars %>% bind_rows(age_mean)


#### ethnicity -----------------------------------------------------------------
## full ethnicity coding
ethnicity <- pair_cc_analytic %>% group_by(ethn_dv_t0) %>% 
  summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Ethnicity") %>% 
  rename("measure"= "ethn_dv_t0") %>% 
  dplyr::select(var, measure, n, est)

## white/non-white
white_non <- pair_cc_analytic %>% group_by(non_white_t0) %>% 
  summarise(n=n()) %>%  
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Ethnicity") %>% 
  rename("measure"= "non_white_t0") %>% 
  dplyr::select(var, measure, n, est) %>% 
  arrange(factor(measure, levels = c("White","Non-white","Missing")))


sample_chars <- sample_chars %>% bind_rows(white_non)

#### Marital status -------------------------------------------------------------
marital <- pair_cc_analytic %>% group_by(marital_status_t0) %>% summarise(n=n()) %>%  
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Marital status") %>% 
  rename("measure"= "marital_status_t0") %>% 
  dplyr::select(var, measure, n, est) %>% 
  arrange(factor(measure, levels = c("married/civil partnership","divorced/separated/widowed","single","missing")))


sample_chars <- sample_chars %>% bind_rows(marital)

#### Educational attainment ----------------------------------------------------
ed_attain <- pair_cc_analytic %>% group_by(hiqual_dv_t0) %>% summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Educational attainment") %>% 
  rename("measure"= "hiqual_dv_t0") %>% 
  dplyr::select(var, measure, n, est) %>% 
  arrange(factor(measure, levels = c("degree",
                                           "other higher degree",
                                           "a-level etc",
                                           "gcse etc",
                                           "other qualification",
                                           "no qualification")))


sample_chars <- sample_chars %>% bind_rows(ed_attain)

#### Region --------------------------------------------------------------------
region <- pair_cc_analytic %>% group_by(gor_dv_t0) %>% summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Region") %>% 
  rename("measure"= "gor_dv_t0") %>% 
  dplyr::select(var, measure, n, est)

sample_chars <- sample_chars %>% bind_rows(region)

#####----------------------------------------------------------------------#####
#####               Employment and income characteristics                  #####
#####----------------------------------------------------------------------#####

#### Current job: SIC-2007 -----------------------------------------------------
sic2007 <- pair_cc_analytic %>% group_by(sic2007_section_lab_t0) %>% summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="SIC-2007") %>% 
  rename("measure"= "sic2007_section_lab_t0") %>% 
  dplyr::select(var, measure, n, est)

sample_chars <- sample_chars %>% bind_rows(sic2007)

#### Current job: SOC-2000 -----------------------------------------------------
soc2000 <- pair_cc_analytic %>% group_by(soc2000_major_group_title_t0) %>% summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="SOC-2000") %>% 
  rename("measure"= "soc2000_major_group_title_t0") %>% 
  dplyr::select(var, measure, n, est)

sample_chars <- sample_chars %>% bind_rows(soc2000)

#### permanent or temporary ----------------------------------------------------
perm_emp <- pair_cc_analytic %>% group_by(emp_contract_t0) %>% summarise(n=n()) %>% 
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Employment contract") %>% 
  rename("measure"= "emp_contract_t0") %>% 
  dplyr::select(var, measure, n, est) %>% 
  arrange(factor(measure, levels = c("fixed-term",
                                           "permanent",
                                           "unemployed/not in employment")))



sample_chars <- sample_chars %>% bind_rows(perm_emp)


#### Employment spells since last interview ------------------------------------

#emp_broken <- pair_cc_analytic %>% group_by(broken_emp_t0) %>% 
#  summarise(n=n()) %>%  
#  mutate(est = n/sum(n)*100) %>% 
#  mutate(var="Broken employment") %>% 
#  rename("measure"= "broken_emp_t0") %>% 
#  dplyr::select(var, measure, n, est)  %>% 
#  arrange(factor(measure, levels = c("Unbroken employment",
#                                           "Broken employment",
#                                           "No employment spells")))


#sample_chars <- sample_chars %>% bind_rows(emp_broken)


#### Multiple jobs -------------------------------------------------------------
### has a 2nd job ----
emp_2nd <- pair_cc_analytic %>% group_by(j2has_dv_t0) %>% summarise(n=n()) %>%
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Multiple jobs") %>% 
  rename("measure"= "j2has_dv_t0") %>% 
  dplyr::select(var, measure, n, est)  %>% 
  arrange(factor(measure, levels = c("no",
                                           "yes")))



sample_chars <- sample_chars %>% bind_rows(emp_2nd)

#### income --------------------------------------------------------------------
## total net personal income (check what used in COVID modelling)
# check monthly?
pair_cc_analytic$fimnnet_dv_t0 <- as.numeric(as.character(pair_cc_analytic$fimnnet_dv_t0))


inc_quantile <- pair_cc_analytic %>%  
  summarise(enframe(quantile(fimnnet_dv_t0, c(0.25, 0.5, 0.75)), "measure", "est")) %>%
  #  mutate(measure=factor(measure)) %>% 
  mutate(measure = ifelse(measure=="25%","25% quantile", 
                          ifelse(measure=="50%","Median","75% quantile"))) %>% 
  mutate(var="Monthly net income (£)",
         n=NA) %>% 
  dplyr::select(var, measure, n, est)

sample_chars <- sample_chars %>% bind_rows(inc_quantile)


#####----------------------------------------------------------------------#####
#####                        Health characteristics                        #####
#####----------------------------------------------------------------------#####

#### long-standing illness or impairment ---------------------------------------
ltc <- pair_cc_analytic %>% group_by(health_t0) %>% summarise(n=n()) %>%  
  mutate(pc = n/sum(n)*100) %>%
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Long-standing illness or impairment") %>% 
  rename("measure"= "health_t0") %>% 
  dplyr::select(var, measure, n, est)  %>% 
  arrange(factor(measure, levels = c("no",
                                     "yes")))

sample_chars <- sample_chars %>% bind_rows(ltc)

#### self-rated health ---------------------------------------------------------

pair_cc_analytic$srh_dv_t0 <- factor(pair_cc_analytic$srh_dv_t0, 
                           levels = c("excellent", "very good", "good",
                                      "fair", "poor"))

srh <- pair_cc_analytic %>% 
  group_by(srh_dv_t0) %>% summarise(n=n()) %>%  
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Self-rated health") %>% 
  rename("measure"= "srh_dv_t0") %>% 
  dplyr::select(var, measure, n, est)

sample_chars <- sample_chars %>% bind_rows(srh)


#### GHQ-12 --------------------------------------------------------------------

ghq4 <- pair_cc_analytic %>% group_by(ghq_case4_t0) %>% summarise(n=n()) %>%  
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="GHQ12 score") %>% 
  rename("measure"= "ghq_case4_t0") %>% 
  dplyr::select(var, measure, n, est) %>% 
  arrange(factor(measure, levels = c("0-3",
                                           "4 or more")))


sample_chars <- sample_chars %>% bind_rows(ghq4)

#### SF-12 physical component summary -------------------------------------------- 

# convert to numeric (to character stirng first to actual value is retained)
pair_cc_analytic$sf12pcs_dv_t0 <- as.character(pair_cc_analytic$sf12pcs_dv_t0)
pair_cc_analytic$sf12pcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12pcs_dv_t0)

# distribution measures
mean(pair_cc_analytic$sf12pcs_dv_t0, na.rm = TRUE)
median(pair_cc_analytic$sf12pcs_dv_t0, na.rm = TRUE)
min(pair_cc_analytic$sf12pcs_dv_t0, na.rm = TRUE)
max(pair_cc_analytic$sf12pcs_dv_t0, na.rm = TRUE)
sf12pcs_quantile <- quantile(pair_cc_analytic$sf12pcs_dv_t0, na.rm = TRUE)

# add to df

#### SF-12 mental component summary -------------------------------------------- 

# convert to numeric (to character stirng first to actual value is retained)
pair_cc_analytic$sf12mcs_dv_t0 <- as.character(pair_cc_analytic$sf12mcs_dv_t0)
pair_cc_analytic$sf12mcs_dv_t0 <- as.numeric(pair_cc_analytic$sf12mcs_dv_t0)

# distribution measures
mean(pair_cc_analytic$sf12mcs_dv_t0, na.rm = TRUE)
median(pair_cc_analytic$sf12mcs_dv_t0, na.rm = TRUE)
min(pair_cc_analytic$sf12mcs_dv_t0, na.rm = TRUE)
max(pair_cc_analytic$sf12mcs_dv_t0, na.rm = TRUE)
sf12mcs_quantile <- quantile(pair_cc_analytic$sf12mcs_dv_t0, na.rm = TRUE)

# add to df


#####----------------------------------------------------------------------#####
#####                          Save sample chars data                      #####
#####----------------------------------------------------------------------#####

## as dataframe
write_rds(sample_chars, "./working_data/sample_chars.rds")



#### job loss between t0-t1

emp_status_t1 <- pair_cc_analytic %>%
  group_by(exposure1)  %>%  
  summarise(n=n()) %>%  
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Employment status t1") %>% 
  rename("measure"= "exposure1") %>% 
  dplyr::select(var, measure, n, est)

# if unemployed at t0 or if unemp spell between t0-1
lost_job <- pair_cc_analytic %>% 
  group_by(exposure2) %>% 
  summarise(n=n()) %>%  
  mutate(est = n/sum(n)*100) %>% 
  mutate(var="Job loss between t0 and t1") %>% 
  rename("measure"= "exposure2") %>% 
  dplyr::select(var, measure, n, est)

################################################################################
#####        Table 1: participant characteristics by exposure group        #####
################################################################################

pair_cc_analytic <- droplevels(pair_cc_analytic)

pair_cc_analytic <- pair_cc_analytic %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(.x%in%c("missing","Missing"),NA,.x))) 

cov_vector <- c("sex_dv_t0", 
                "age_dv_t0",  
                "non_white_t0", 
                "marital_status_t0",
                "hiqual_dv_t0", 
                "gor_dv_t0",
                "sic2007_section_lab_t0",
                "soc2000_major_group_title_t0",
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

## unemployed at T1
table_one <- CreateTableOne(vars = cov_vector, strata = "exposure1",
                            data = pair_cc_analytic,
                            test = FALSE)

table_one_sav <- print(table_one, showAllLevels = TRUE, smd = TRUE, 
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

