################################################################################

# Precarious employment and health - Understanding Society
# 1-03 - Clean and recode master raw data for causal analysis
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) .
# (b) .  
# (c) .
# Data output: cleaned master raw dataframe

#### To be added:

################################################################################

## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(readxl) # for reading excel files
library(foreign) # for reading SPSS files
library(tidyverse) # all kinds of stuff 
library(janitor) # cleaning up

################################################################################
#####                               load functions                         #####
################################################################################

## function for creating occupation and industry group vars
source("./functions/soc_sic_prep_function.R")


################################################################################
#####                               load data                              #####
################################################################################

master_raw1 <- readRDS("./raw_data/master_raw1.rds")




################################################################################
#####                               identifiers                           ######
################################################################################

master_raw1 <- master_raw1 %>% 
  mutate(retchk_flag = ifelse(retchk %in% c("yes, considers is retired",
                                       "Yes, considers is retired"),1,0))

################################################################################
#####                                covariates                            #####
################################################################################

### recode inconsistent sex as missing
master_raw1 <- master_raw1 %>% 
  mutate(sex_dv = as.character(sex_dv)) %>% 
  mutate(sex_dv = ifelse(sex_dv=="inconsistent","missing",sex_dv))

### recode age as numeric
master_raw1 <- master_raw1 %>% 
  mutate(age_dv=as.character(age_dv)) %>% 
  mutate(age_dv=as.numeric(age_dv))


### recode ethnicity as white/non-white
master_raw1 <- master_raw1 %>% 
  mutate(non_white = ifelse(ethn_dv=="british/english/scottish/welsh/northern irish", "White",
                            ifelse(ethn_dv=="irish", "White",
                                   ifelse(ethn_dv=="any other white background", "White",
                                          ifelse(ethn_dv=="missing", "Missing",
                                                 "Non-white")))))


### marital status
master_raw1$marstat <- as.character(master_raw1$marstat)
master_raw1 <- master_raw1 %>% 
  mutate(marital_status = ifelse(marstat %in% c("missing",                                                     
                                                "inapplicable",                                                
                                                "proxy",                                                       
                                                "refusal",                                                     
                                                "don't know",                                                  
                                                "Only available for IEMB",                                     
                                                "Not available for IEMB"),
                                 "missing", 
                           ifelse(marstat %in% c("single and never married or never in a legally recognised ci",
                                                "single & never married or in legally recog'd civil p'ship",
                                                "single, nvr marr/civ p",                                      
                                                "Single, nvr marr/civ p"),
                                 "single",
                           ifelse(marstat %in% c("married",                                                     
                                                "a civil partner in a legally recognised civil partnership",   
                                                "civil partner (legal)",                                       
                                                "Married",                                                     
                                                "Civil Partner (legal)",
                                                "in a registered same-sex civil partnership"),
                                  "married/civil partnership",
                           ifelse(marstat %in% c("separated but legally married",                               
                                                 "divorced",                                                    
                                                 "widowed", 
                                                 "spontaneous: separated from civil partner",                   
                                                 "spontaneous: a former civil partner, the civil partnership l",
                                                 "spontaneous: a surviving civil partner (your partner having", 
                                                 "separated legally marr",                                      
                                                 "sep from civil partner",                                      
                                                 "a former civil partner",                                      
                                                 "surviving civil partner",                                     
                                                 "Separated legally marr",                                      
                                                 "Divorced",                                                    
                                                 "Widowed",                                                     
                                                 "Sep from Civil Partner",                                      
                                                 "A former Civil Partner",                                      
                                                 "Surviving Civil Partner",
                                                 "an ex-civil partner,civil p'ship legally dissolved",
                                                 "surviving civil partner (partner died)",
                                  "separated from civil partner"),
                                  "divorced/separated/widowed",
                                  "check")))))
                                      

### educational attainment

master_raw1$hiqual_dv <- as.character(master_raw1$hiqual_dv)

master_raw1 <- master_raw1 %>% 
  mutate(hiqual_dv = ifelse(hiqual_dv=="Other higher","Other higher degree",
                  ifelse(hiqual_dv=="Other qual","Other qualification",
                         ifelse(hiqual_dv=="No qual","No qualification",
                                ifelse(hiqual_dv=="A level etc","A-level etc",
                                       ifelse(hiqual_dv=="inapplicable","missing",
                                       hiqual_dv))))))

  
### region
master_raw1 <- master_raw1 %>% 
  mutate(gor_dv = tolower(gor_dv))
           
### perceived job security
master_raw1 <- master_raw1 %>% 
  mutate(jbsec_dv = ifelse(jbsec %in% c("missing",                 
                         "inapplicable",            
                         "proxy",                   
                         "refusal",                
                         "don't know",            
                         "Only available for IEMB", 
                         "Not available for IEMB"),
         "missing",
         ifelse(jbsec %in% c("very likely",
                             "Very likely"),
                "very likely",
                ifelse(jbsec %in% c("Likely",                 
                                    "likely"),
                       "likely",
                       ifelse(
                         jbsec %in% c("unlikely",               
                                      "Unlikely"),
                         "unlikely",
                         ifelse(jbsec %in% c("very unlikely?",    
                                             "Very unlikely?"),
                                "very unlikely",
                                "missing"))))))
  


### employment

master_raw1$employ <- as.character(tolower(master_raw1$employ))

master_raw1 <- master_raw1 %>% 
  mutate(employ = ifelse(employ %in% "employer moved job to another workplace",
                         "yes",
                         ifelse(employ %in% "got a different job with the same employer which meant movin",
                                "no",employ)))


 
### long-standing condition
master_raw1 <- master_raw1 %>% 
  mutate(health = ifelse(health %in% c("yes","Yes"),"yes",
                  ifelse(health %in% c("no","No"),"no",
                  ifelse(health %in% c("missing","refusal",
                                       "don't know","inapplicable"),"missing",
                         "check"))))

### industry (SIC-2007)
master_raw1 <- sic2007f()

## create "other" category for smaller industries (n<10k)
master_raw1 <- master_raw1 %>%
  mutate(sic2007_section_lab = 
           ifelse(sic2007_section_lab %in% 
                    c("Activities of extraterritorial organisations and bodies",
                      "Activities of households as employers; undifferentiated goods-and services-producing activities of households for own use",
                      "Agriculture, forestry and fishing",
                      "Arts, entertainment and recreation",
                      "Electricity, gas, steam and air conditioning supply",
                      "Financial and insurance activities",
                      "Information and communication",
                      "Mining and quarrying",
                      "Other service activities",
                      "Real estate activities",
                      "Water supply; sewerage, waste management and remediation activities"),
                  "Other industry", sic2007_section_lab))
  
master_raw1 %>% 
  group_by(sic2007_section_lab) %>% 
  summarise(n=n()) %>% 
  print(n=25)

### occupation (SOC-2000)
master_raw1 <- soc2000f()

### normal weekly hours
## convert to character rather than factor
master_raw1$jbhrs <- as.character(master_raw1$jbhrs)

## change missing categories to NA then convert to numeric
master_raw1 <- master_raw1 %>% 
  mutate(jbhrs = ifelse(jbhrs %in% c("missing", "inapplicable", "proxy", 
                                     "refusal", "don't know"),NA, jbhrs))
master_raw1$jbhrs <- as.numeric(master_raw1$jbhrs)

### income
master_raw1$fimnnet_dv <- as.character(master_raw1$fimnnet_dv)
master_raw1$fimnnet_dv <- as.numeric(master_raw1$fimnnet_dv)


## check <0s
master_raw1 %>% mutate(inc_flag = ifelse(fimnnet_dv<0,1,0)) %>% 
  group_by(inc_flag) %>% 
  summarise(n=n())
# can be <0 due to self-employ,ent losses reported
# so shouldn't see any/many in analytic sample

## calculate median income and relative poverty for each wave

inc_temp <- master_raw1 %>% 
  group_by(wv_n) %>% 
  summarise(med_inc = median(fimnnet_dv, na.rm=TRUE)) %>% 
  mutate(med_inc_60 = 0.6*med_inc) # calculate 60% median income (relative poverty)

# join back onto master df
master_raw1 <- master_raw1 %>% left_join(inc_temp) %>% 
  mutate(below_med_inc = ifelse(fimnnet_dv<med_inc,"below median income","median income of above"),
         rel_pov = ifelse(fimnnet_dv<med_inc_60,"relative poverty","not in relative poverty")) %>% 
  dplyr::select(-c(med_inc, med_inc_60))

table(master_raw1$below_med_inc)
table(master_raw1$rel_pov)


################################################################################
#####                            exposure variables                        #####
################################################################################

#### recode employment status variables to create an employment contract variable
master_raw1 <- master_raw1 %>% 
  mutate(emp_contract = ifelse(jbterm1 %in% c("not permanent job", 
                                              "or is there some way that it is not permanent?",
                                              "Or is there some way that it is not permanent?"),
                               "fixed-term",
                               ifelse(jbterm1 %in% c("a permanent job",
                                                     "A permanent job"),
                                      "permanent",
                                      ifelse(jbstat %in% c("Unemployed",
                                                           "unemployed",
                                                           "on maternity leave",
                                                           "On maternity leave",
                                                           "Family care or home",
                                                           "full-time student",
                                                           "Full-time student",
                                                           "LT sick or disabled",
                                                           "Govt training scheme", # or fixed-term?
                                                           "On apprenticeship",    # or fixed-term?
                                                           "Unpaid, family business",
                                                           "doing something else",
                                                           "Doing something else"), 
                                                    "unemployed/not in employment",
                                                    "missing"))))

#### recode employment spells vars to create a broken employment variable
### employment spells -----
# change to character var
master_raw1$nmpsp_dv <- as.character(master_raw1$nmpsp_dv)
master_raw1 <- master_raw1 %>% 
  mutate(nmpsp_dv=ifelse(nmpsp_dv %in% c( "missing", "inapplicable", "proxy", "refusal", "don't know"),"missing",
                       ifelse(nmpsp_dv=="none","0",nmpsp_dv))) 

master_raw1 <- master_raw1 %>% 
  mutate(emp_spells_bin = ifelse(nmpsp_dv == "missing", "missing",
                          ifelse(nmpsp_dv == "0","no",
                                 "yes")))

# checks
table(master_raw1$nmpsp_dv,master_raw1$emp_spells_bin)
table(master_raw1$emp_spells_bin,master_raw1$wv_n)

sum(master_raw1$emp_spells_bin=="missing")

### non-employment spells -----
# change to character var
master_raw1$nnmpsp_dv <- as.character(master_raw1$nnmpsp_dv)

master_raw1 <- master_raw1 %>% 
  mutate(nnmpsp_dv=ifelse(nnmpsp_dv %in% c( "missing", "inapplicable", "proxy", "refusal", "don't know"),"missing",
                          ifelse(nnmpsp_dv=="none","0",nnmpsp_dv)))

master_raw1 <- master_raw1 %>% 
  mutate(non_emp_spells_bin = ifelse(nnmpsp_dv == "missing", "missing",
                                 ifelse(nnmpsp_dv == "0","no",
                                        "yes")))

# checks
table(master_raw1$nnmpsp_dv,master_raw1$non_emp_spells_bin)
table(master_raw1$non_emp_spells_bin,master_raw1$wv_n)

sum(master_raw1$non_emp_spells_bin=="missing")


### unemployment spells -----
# change to character var
master_raw1$nunmpsp_dv <- as.character(master_raw1$nunmpsp_dv)

master_raw1 <- master_raw1 %>% 
  mutate(nunmpsp_dv=ifelse(nunmpsp_dv %in% c( "missing", "inapplicable", "proxy", "refusal", "don't know"),"missing",
                          ifelse(nunmpsp_dv=="none","0",nunmpsp_dv)))

master_raw1 <- master_raw1 %>% 
  mutate(unemp_spells_bin = ifelse(nunmpsp_dv == "missing", "missing",
                                     ifelse(nunmpsp_dv == "0","no",
                                            "yes")))

# checks
table(master_raw1$nunmpsp_dv,master_raw1$unemp_spells_bin)
table(master_raw1$unemp_spells_bin,master_raw1$wv_n)

sum(master_raw1$unemp_spells_bin=="missing")

### create broken employment variable ------
master_raw1 <- master_raw1 %>% 
  mutate(broken_emp = ifelse(emp_spells_bin=="no","No employment spells", 
                             ifelse(emp_spells_bin=="yes" & 
                                      non_emp_spells_bin=="no" & 
                                      unemp_spells_bin=="no",
                                    "Unbroken employment",
                                    ifelse(emp_spells_bin=="yes" & 
                                             (non_emp_spells_bin=="yes" |
                                                unemp_spells_bin=="yes"),
                                           "Broken employment","missing"))))

#### recode multiple jobs var for analysis

## collapse missing categories
master_raw1$j2has <- as.character(master_raw1$j2has)


master_raw1 <- master_raw1 %>% 
  mutate(j2has_dv = ifelse(j2has %in% c("proxy","missing","refusal","don't know"),
                           "missing",
                           ifelse(j2has %in% c("Yes","yes"),
                                  "yes",
                                  "no")))

master_raw1 <- master_raw1 %>% 
  mutate(j2has_dv2 = ifelse(j2has_dv == "yes" & 
                              emp_contract %in% c("fixed-term","permanent"), 
                            "multiple jobs",
                            ifelse(j2has_dv == "yes" & 
                                     emp_contract %in% c("unemployed/not in employment","missing"),
                                   "unemployed/not in employment with additional",
                            ifelse(j2has_dv %in% c("no","missing") & 
                                     emp_contract %in% c("fixed-term","permanent"), 
                                   "one job",
                                   ifelse(j2has_dv %in% c("no","missing") & 
                                            emp_contract=="unemployed/not in employment",
                                   "unemployed/not in employment",
                                   "missing")))))


#### full/part-time employment recode
master_raw1$jbft_dv <- as.character(master_raw1$jbft_dv)

master_raw1 <- master_raw1 %>% 
  mutate(jbft_dv = ifelse(jbft_dv %in% c("proxy","missing","refusal","don't know", "inapplicable"),
                          "missing",
                          jbft_dv))

#### employer size recode
master_raw1$jbsize <- as.character(master_raw1$jbsize)


master_raw1 <- master_raw1 %>% 
  mutate(jbsize = ifelse(jbsize %in% c("proxy","missing","refusal","don't know", "inapplicable"),
                         "missing",
                         ifelse(jbsize == "Don't know but fewer than 25","10 - 24",
                                ifelse(jbsize == "Don't know but 25 or more","25 - 49",
                         jbsize)))) %>% 
  mutate(small_firm = ifelse(jbsize=="missing","missing",
                             ifelse(jbsize %in% c("1 - 2",  "10 - 24", "25 - 49"),
                                    "under 50 employees",
                                    "over 50 employees")))



################################################################################
#####                            health outcomes                           #####
################################################################################

#### recode self-rated health for analysis

master_raw1 <- master_raw1 %>% 
  # recode self-rated health variables into one
  mutate(sf1 = as.character(sf1),
         scsf1 = as.character(scsf1)) %>% 
  mutate(srh_dv = ifelse(sf1=="inapplicable",scsf1, sf1)) %>% 
  mutate(srh_dv = ifelse(srh_dv %in% c("excellent", "Excellent"),
                         "excellent",
                  ifelse(srh_dv %in% c("very good", "Very good"),
                         "very good",
                  ifelse(srh_dv %in% c( "good","Good"),
                         "good",
                  ifelse(srh_dv %in% c("fair", "Fair"),
                         "fair",
                  ifelse(srh_dv %in% c("poor", "Poor", "or Poor?"),
                         "poor",
                  ifelse(srh_dv %in% c("don't know", "inapplicable",
                                       "missing", "proxy",  "refusal"),
                         "missing",
                         "check"
                         )))))))
           
## binary version for outcome analysis
master_raw1 <- master_raw1 %>% 
  mutate(srh_bin = ifelse(srh_dv %in% c("excellent","very good"),
                          "excellent/very good",
                   ifelse(srh_dv %in% c("good","fair","poor"),
                          "good/fair/poor",
                          "missing")))


#### recode GHQ-12 caseness for analysis
## calculate caseness for main analysis (cut point = 4)
master_raw1 <- master_raw1 %>% 
  mutate(ghq_case4 = ifelse(grepl("0",as.character(scghq2_dv)),0,
                          ifelse(grepl("1",as.character(scghq2_dv)),0,
                                 ifelse(grepl("2",as.character(scghq2_dv)),0,
                                        ifelse(grepl("3",as.character(scghq2_dv)),0,
                                               ifelse(grepl("4",as.character(scghq2_dv)),1,
                                                      ifelse(grepl("5",as.character(scghq2_dv)),1,
                                                             ifelse(grepl("6",as.character(scghq2_dv)),1,
                                                                    ifelse(grepl("7",as.character(scghq2_dv)),1,
                                                                           ifelse(grepl("8",as.character(scghq2_dv)),1,
                                                                                  ifelse(grepl("9",as.character(scghq2_dv)),1,
                                                                                         ifelse(grepl("10",as.character(scghq2_dv)),1,
                                                                                                ifelse(grepl("11",as.character(scghq2_dv)),1,
                                                                                                       ifelse(grepl("12",as.character(scghq2_dv)),1,
                                                                                                              as.character(scghq2_dv))))))))))))))) 

master_raw1 <- master_raw1 %>% 
  mutate(ghq_case4 = ifelse(ghq_case4=="0","0-3",
                            ifelse(ghq_case4=="1","4 or more",
                                   ifelse(ghq_case4 %in% c("inapplicable", 
                                                           "missing", 
                                                           "proxy"),
                                          "missing",
                                          ghq_case4))))

## calculate caseness for sensitivity analysis (cut point = 3)
master_raw1 <- master_raw1 %>% mutate(ghq_case3 = ifelse(grepl("0",as.character(scghq2_dv)),0,
                                                         ifelse(grepl("1",as.character(scghq2_dv)),0,
                                                                ifelse(grepl("2",as.character(scghq2_dv)),0,
                                                                       ifelse(grepl("3",as.character(scghq2_dv)),1,
                                                                              ifelse(grepl("4",as.character(scghq2_dv)),1,
                                                                                     ifelse(grepl("5",as.character(scghq2_dv)),1,
                                                                                            ifelse(grepl("6",as.character(scghq2_dv)),1,
                                                                                                   ifelse(grepl("7",as.character(scghq2_dv)),1,
                                                                                                          ifelse(grepl("8",as.character(scghq2_dv)),1,
                                                                                                                 ifelse(grepl("9",as.character(scghq2_dv)),1,
                                                                                                                        ifelse(grepl("10",as.character(scghq2_dv)),1,
                                                                                                                               ifelse(grepl("11",as.character(scghq2_dv)),1,
                                                                                                                                      ifelse(grepl("12",as.character(scghq2_dv)),1,
                                                                                                                                             as.character(scghq2_dv))))))))))))))) 

master_raw1 <- master_raw1 %>% 
  mutate(ghq_case3 = ifelse(ghq_case3=="0","0-2",
                            ifelse(ghq_case3=="1","3 or more",
                                   ifelse(ghq_case3 %in% c("inapplicable", 
                                                           "missing", 
                                                           "proxy"),
                                          "missing",
                                          ghq_case3))))


#### sf-12 PCS
master_raw1 <- master_raw1 %>% 
  mutate(sf12pcs_dv = as.character(sf12pcs_dv)) %>% 
  mutate(sf12pcs_dv = ifelse(sf12pcs_dv %in%c("inapplicable", 
                                              "missing", 
                                              "proxy",
                                              "(Other)"),
                             "missing", sf12pcs_dv))

#### sf-12 MCS
master_raw1 <- master_raw1 %>% 
  mutate(sf12mcs_dv = as.character(sf12mcs_dv)) %>% 
  mutate(sf12mcs_dv = ifelse(sf12mcs_dv %in%c("inapplicable", 
                                              "missing", 
                                              "proxy",
                                              "(Other)"),
                             "missing", sf12mcs_dv))


################################################################################
#####                         save cleaned dataframe                       #####
################################################################################

## write as RDS
write_rds(master_raw1, "./working_data/master_raw1_clean.rds")

