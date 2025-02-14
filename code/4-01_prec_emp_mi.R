################################################################################

# Precarious employment and health - Understanding Society
# 4-01 - Create multiple imputed analytic sample 
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Checks item missingness across survey waves 
# (b) create final MI df 
# (c) produce unweighted/matched descriptives and outcomes analysis


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
library(tableone)

`%notin%` <- Negate(`%in%`)

################################################################################
#####                         load and prepare data                        #####
################################################################################


####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")


####load eligible cases --------------------------------------------------------
pair_eligible <- readRDS("./working_data/pair_eligible.rds") %>% 
    dplyr::select(pidp, exposure1, all_of(c(cov_vector, cov_vector2, 
                           outcome_vector2)))



# check number of missing cases by variable
sapply(pair_eligible, function(x) sum(is.na(x)))

# check number of missing cases by row

temp <- pair_eligible %>% 
  dplyr::select(# person identifier
    pidp,
    # outcomes vars
    sf12pcs_dv_t1,
    sf12mcs_dv_t1, 
    srh_bin_t1,
    ghq_case4_t1,
    # exposure 
    exposure1, 
    # t0 covariates
    all_of(cov_vector),
    # t1 covariates
    age_dv_t1,
    marital_status_t1,
    health_t1
  )
temp$na_count <- rowSums(is.na(temp))

table(temp$na_count)


## calculate propotion of rows with missing data - use for number of imputations
df_rows <- nrow(pair_eligible)

cc_rows <- nrow(na.omit(pair_eligible))

1-cc_rows/df_rows

rm(temp)

################################################################################
#####                             create NAs df                            #####
################################################################################

#missing_vector <- c("missing", "Missing",
#                    "inapplicable", "proxy",
#                    "refusal", 
#                    "Only available for IEMB", 
#                    "Not available for IEMB",
#                    "don't know")
#

# check
sapply(pair_eligible, function(x) sum(is.na(x)))

### convert string missing terms into NAs
pair_eligible <- pair_eligible %>% mutate(across(.cols = everything(),
                                        .fns = ~ifelse(.%in% c("missing", "Missing",
                                                 "inapplicable", "proxy",
                                                 "refusal", 
                                                 "Only available for IEMB", 
                                                 "Not available for IEMB",
                                                 "don't know"),NA,.x)))

# check again
sapply(pair_eligible, function(x) sum(is.na(x)))

### create NA df
pair_eligible_na <- pair_eligible %>% 
  #  mutate(across(everything(), as.character)) %>% 
  mutate(across(.cols = everything(), 
                .fns = ~ifelse(is.na(.x),1,0)))



################################################################################
#####                       item missing descriptives                      #####
################################################################################

#### calculate total number of missing cases by variable -----------------------
n_row <- nrow(pair_eligible_na)

na_by_var <- pair_eligible_na %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(.x))) %>% 
  pivot_longer(cols=everything(), names_to = "variable", values_to = "n_NA") %>% 
  mutate(pc_NA = n_NA/n_row*100) 

#### calculate total number of missing items per case --------------------------
na_by_case <- pair_eligible
na_by_case$n_NA <- apply(X = is.na(pair_eligible), MARGIN = 1, FUN = sum)
na_by_case <- na_by_case %>% dplyr::select(pidp, n_NA)

summary(na_by_case$n_NA)

hist(na_by_case$n_NA)

#### visualise missingness -----------------------------------------------------
### create a dataframe with only vars that have NAs
incomplete_vars <- pair_eligible[,colSums(is.na(pair_eligible))>0]

### create missingness pattern plot
tiff("./output/mi/mi_descriptives/md_pattern.tiff",
     width = 6000, height = 32000, res = 300)
plot_pattern( incomplete_vars, vrb="all", square=TRUE, 
              rotate=TRUE, cluster=NULL, npat=NULL, caption=FALSE ) +
  theme(text=element_text(size=14))
dev.off()

### create missingness correlation matrix
cormat <- round(cor(pair_eligible_na),2)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
lower_tri <- get_lower_tri(cormat)

# vector of vars with no NAs
no_NAs <- c("exposure1",
            "exposure2",
            "retchk_flag_t1",
            "jbstat_t0",
            "age_dv_t1",
            "rel_pov_t0",
            "pidp",
            "age_dv_t0",
            "nunmpsp_dv_t0",
            "dep_child_bin_t0")

# melt into format for ggplot
cor_melt <- reshape2::melt(lower_tri) %>% 
  replace(is.na(.), 0) %>% 
# remove cases for no missing vars
    filter(!(Var1%in%no_NAs) & !(Var2%in%no_NAs))
  
## make the plot
ggplot(data = cor_melt, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, 
                                   size = 8, hjust = 1))+
  theme(axis.text.y = element_text(size = 8))+
  coord_fixed()

#####   create interaction terms between outcome vars and sub-group vars   #####

### convert outcomes to numeric
## SF-12 PCS
pair_eligible$sf12pcs_dv_t0 <- as.character(pair_eligible$sf12pcs_dv_t0)
pair_eligible$sf12mcs_dv_t0 <- as.character(pair_eligible$sf12mcs_dv_t0)
pair_eligible$sf12pcs_dv_t1 <- as.character(pair_eligible$sf12pcs_dv_t1)
pair_eligible$sf12mcs_dv_t1 <- as.character(pair_eligible$sf12mcs_dv_t1)

# t1 score - mean for interaction term
#pair_eligible <- pair_eligible %>% 
#  mutate(sf12pcs_dv_t1i = sf12pcs_dv_t1-mean(sf12pcs_dv_t1))

## SF-12 MCS
pair_eligible$sf12pcs_dv_t0 <- as.numeric(pair_eligible$sf12pcs_dv_t0)
pair_eligible$sf12mcs_dv_t0 <- as.numeric(pair_eligible$sf12mcs_dv_t0)
pair_eligible$sf12pcs_dv_t1 <- as.numeric(pair_eligible$sf12pcs_dv_t1)
pair_eligible$sf12mcs_dv_t1 <- as.numeric(pair_eligible$sf12mcs_dv_t1)

# t1 score - mean for interaction term
#pair_eligible <- pair_eligible %>% 
#  mutate(sf12mcs_dv_t1i = sf12mcs_dv_t1-mean(sf12mcs_dv_t1))

## SRH
pair_eligible <- pair_eligible %>% 
  mutate(srh_bin2 = ifelse(srh_bin_t1 == "excellent/very good",0,
                          ifelse(srh_bin_t1 == "good/fair/poor",1,-99)))

## GHQ-12
pair_eligible <- pair_eligible %>% 
  mutate(ghq_bin = ifelse(ghq_case4_t1 == "0-3",0,
                          ifelse(ghq_case4_t1 == "4 or more",1,-99)))

### covert sub-group vars to bin vars
pair_eligible <- pair_eligible %>% 
  ## sex
  mutate(sex_bin = ifelse(sex_dv_t0=="Female",0,
                          ifelse(sex_dv_t0=="Male",1,-99))) %>% 
  ## age - for sub-group analysis rather than MI model
  mutate(age_bin = ifelse(age_dv_t0<median(age_dv_t0),0,
                          ifelse(age_dv_t0>=median(age_dv_t0),1,-99))) %>% 
  ## relative poverty
  mutate(rel_pov_bin = ifelse(rel_pov_t0=="not in relative poverty",0,
                              ifelse(rel_pov_t0=="relative poverty",1,-99)))

#### sex -----------------------------------------------------------------------

pair_eligible <- pair_eligible %>% 
  mutate(sex_pcs = sex_bin*(sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)), ### SF-12 PCS
         sex_mcs = sex_bin*(sf12mcs_dv_t1-mean(sf12mcs_dv_t1, na.rm = TRUE)),### SF-12 MCS
         sex_srh = sex_bin*srh_bin2, ### SRH
         sex_ghq = sex_bin*ghq_bin) ### GHQ-12

#### age -----------------------------------------------------------------------

pair_eligible <- pair_eligible %>% 
  mutate(age_pcs = age_bin*(sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)), ### SF-12 PCS
         age_mcs = age_bin*(sf12mcs_dv_t1-mean(sf12mcs_dv_t1, na.rm = TRUE)),### SF-12 MCS
         age_srh = age_bin*srh_bin2, ### SRH
         age_ghq = age_bin*ghq_bin) ### GHQ-12

#### relative poverty ----------------------------------------------------------

pair_eligible <- pair_eligible %>% 
  mutate(rel_pov_pcs = rel_pov_bin*(sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)), ### SF-12 PCS
         rel_pov_mcs = rel_pov_bin*(sf12mcs_dv_t1-mean(sf12mcs_dv_t1, na.rm = TRUE)),### SF-12 MCS
         rel_pov_srh = rel_pov_bin*srh_bin2, ### SRH
         rel_pov_ghq = rel_pov_bin*ghq_bin) ### GHQ-12

#### interaction variable vector -----------------------------------------------

interactions_vars <- c("srh_bin2","ghq_bin",
                       "sex_pcs","sex_mcs","sex_srh","sex_ghq",
                       "age_pcs","age_mcs","age_srh","age_ghq",
                       "rel_pov_pcs","rel_pov_mcs","rel_pov_srh","rel_pov_ghq")

################################################################################
#####                Multiple imputation for one variable                  #####
################################################################################

#### create sub-set of variables and cases -------------------------------------
mi_subset1 <- pair_eligible %>% 
  dplyr::select(sf12pcs_dv_t1, sex_dv_t0, age_dv_t0, rel_pov_bin,
                rel_pov_pcs, exposure1)

# check var structure is OK for glm
str(mi_subset1)

# change outcome var to numeric
mi_subset1$sf12pcs_dv_t1 <- as.numeric(mi_subset1$sf12pcs_dv_t1)
# change exposure to factor so it is ordered correctly for model
mi_subset1$exposure1 <- factor(mi_subset1$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))
# change sex to factor
mi_subset1$sex_dv_t0 <- factor(mi_subset1$sex_dv_t0,
                               levels = c("Female",
                                          "Male"))

# take a random sample for testing code (if it's running slow)
mi_subset1 <- sample_n(mi_subset1, 3000)


#### create imputations --------------------------------------------------------

## set seed
set.seed(7341)

## imputations
imps <- mice(mi_subset1, m=15, maxit=10)
summary(imps)

## log reg on SF-12 PCS with one covariate
fit <- with(imps, exp=glm(sf12pcs_dv_t1~sex_dv_t0, family="gaussian"))
summary(pool(fit), conf.int = TRUE)

### convert data from mids format to standard df
temp <- complete(imps, "long", include = FALSE)

## sum number of NAs by var
temp %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x))))

## compare with non-imputed data
mi_subset1  %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x))))

################################################################################
#####             Test passive imputation for sub-group analysis           #####
################################################################################

#### SF-12(PCS) ---------------------------------------------------------------

## imputations
imps_passive_pcs <- mice(mi_subset1,maxit=0)
summary(imps_passive_pcs)

# specify method for each incomplete variable in data
# (numeric, binary, unordered, ordered)
# ("norm", "logreg", "polyreg", "polr")
meth_pcs <- imps_passive_pcs$meth
meth_pcs["sex_dv_t0"] <- "logreg"
meth_pcs["age_dv_t0"] <- "norm"
meth_pcs["exposure1"] <- "logreg"
meth_pcs["rel_pov_bin"] <- "pmm"
meth_pcs["sf12pcs_dv_t1"] <- "norm"
meth_pcs["rel_pov_pcs"] <- "~I(rel_pov_bin*(sf12pcs_dv_t1-mean(sf12pcs_dv_t1)))"

pred_pcs_passive <- imps_passive_pcs$pred
pred_pcs_passive[c("rel_pov_bin","sf12pcs_dv_t1"),"rel_pov_pcs"] <- 0

visit_pcs <- imps_passive_pcs$visitSequence
visit_pcs2 <- c("sex_dv_t0", "age_dv_t0", "exposure1", "rel_pov_bin", 
            "sf12pcs_dv_t1", "rel_pov_pcs")

imps_passive_pcs2 <- mice(mi_subset1, pred = pred_pcs_passive, meth = meth_pcs, maxit = 10, 
                      visitSequence = visit_pcs2)

## log reg on SF-12 PCS with one covariate
fit_pcs <- with(imps_passive_pcs2, exp=glm(sf12pcs_dv_t1~sex_dv_t0, family="gaussian"))
summary(pool(fit_pcs), conf.int = TRUE)

### convert data from mids format to standard df
temp_pcs <- complete(imps_passive_pcs2, "long", include = FALSE)

## sum number of NAs by var
temp_pcs %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x))))

temp_pcs %>% head()

#### check ---------------------------------------------------------------------

check_pcs_df1 <- pair_eligible %>% filter(is.na(rel_pov_bin))

check_pcs_df2 <- pair_eligible  %>% 
  filter(exposure1=="exposed (employed at t1)" & age_dv_t0%in%c(54) & 
           sf12pcs_dv_t1%in%c(47.37))

check_pcs_df3 <- temp_pcs %>% 
  filter(exposure1=="exposed (employed at t1)" & age_dv_t0%in%c(54) & 
           sf12pcs_dv_t1%in%c(47.37))

temp_pcs %>% filter(rel_pov_bin==1) %>% filter(sf12pcs_dv_t1!=rel_pov_pcs)

temp_pcs %>% filter(rel_pov_bin==1) %>% summarise(sum(sf12pcs_dv_t1-rel_pov_pcs))

temp_pcs %>% filter(rel_pov_bin==1) %>% filter(sf12pcs_dv_t1!=rel_pov_pcs) %>% summarise(n=n())

temp_pcs %>% 
  filter(rel_pov_bin==1) %>% 
  ggplot() + 
  geom_point(aes(x = sf12pcs_dv_t1, y = rel_pov_pcs))

#### SRH -----------------------------------------------------------------------

#### create sub-set of variables and cases -------------------------------------
mi_subset_srh <- pair_eligible %>% 
  dplyr::select(srh_bin_t1, sex_dv_t0, age_dv_t0, rel_pov_bin,
                rel_pov_srh, exposure1)

# check var structure is OK for glm
str(mi_subset_srh)

# change outcome var to factor
mi_subset_srh$srh_bin_t1 <- factor(mi_subset_srh$srh_bin_t1,
                                   levels = c("excellent/very good","good/fair/poor"))
# change exposure to factor so it is ordered correctly for model
mi_subset_srh$exposure1 <- factor(mi_subset_srh$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))
# change sex to factor
mi_subset_srh$sex_dv_t0 <- factor(mi_subset_srh$sex_dv_t0,
                               levels = c("Female",
                                          "Male"))

# take a random sample for testing code (if it's running slow)
mi_subset_srh <- sample_n(mi_subset_srh, 3000)

## imputations
imps_passive_srh <- mice(mi_subset_srh,maxit=0)
summary(imps_passive_srh)

# specify method for each incomplete variable in data
# (numeric, binary, unordered, ordered)
# ("norm", "logreg", "polyreg", "polr")
meth_srh <- imps_passive_srh$meth
meth_srh["sex_dv_t0"] <- "logreg"
meth_srh["age_dv_t0"] <- "norm"
meth_srh["rel_pov_bin"] <- "pmm"
meth_srh["exposure1"] <- "logreg"
meth_srh["srh_bin_t1"] <- "pmm"
meth_srh["rel_pov_srh"] <- "~I(rel_pov_bin*srh_bin)"

pred_passive_srh <- imps_passive_srh$pred
pred_passive_srh[c("rel_pov_bin","srh_bin_t1"),"rel_pov_srh"] <- 0

visit_srh <- imps_passive_srh$visitSequence
visit_srh2 <- c("sex_dv_t0", "age_dv_t0", "exposure1", "rel_pov_bin", 
            "sf12pcs_dv_t1", "rel_pov_pcs")

imps_passive_srh <- mice(mi_subset_srh, 
                         pred = pred_passive_srh, 
                         meth = meth_srh, 
                         maxit = 10, 
                      visitSequence = visit_srh2)

## log reg on SF-12 PCS with one covariate
fit_srh <- with(imps_passive_srh, exp=glm(srh_bin_t1~sex_dv_t0, family = binomial(link = logit)))
summary(pool(fit_srh), conf.int = TRUE)

### convert data from mids format to standard df
temp_srh <- complete(imps_passive_srh, "long", include = FALSE)

## sum number of NAs by var
temp_srh %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x))))

temp_srh %>% head()

#### check ---------------------------------------------------------------------

check_srh_df1 <- pair_eligible %>% filter(is.na(rel_pov_bin))

check_srh_df2 <- pair_eligible  %>% 
  filter(exposure1=="exposed (employed at t1)" & age_dv_t0%in%c(54) & 
           sf12pcs_dv_t1%in%c(47.37))

check_srh_df3 <- temp %>% 
  filter(exposure1=="exposed (employed at t1)" & age_dv_t0%in%c(54) & 
           sf12pcs_dv_t1%in%c(47.37))

temp_srh %>% filter(srh_bin_t1==1 & rel_pov_bin==1) %>% filter(rel_pov_srh==1)

temp_srh %>% filter(rel_pov_srh!=1 &rel_pov_srh!=0)

temp_srh %>% 
  filter(rel_pov_bin==1) %>% 
  ggplot() + 
  geom_point(aes(x = srh_bin_t1, y = rel_pov_srh))

################################################################################
#####             Multiple imputation for multiple variables               #####
################################################################################

#### create sub-set of variables and cases -------------------------------------
mi_subset2 <- pair_eligible %>% 
  dplyr::select(# person identifier
                pidp,
                # outcomes vars
                sf12pcs_dv_t1,
                sf12mcs_dv_t1, 
                srh_bin_t1,
                ghq_case4_t1,
                # exposure 
                exposure1, 
                # t0 covariates
                all_of(cov_vector),
                # t1 covariates
                age_dv_t1,
                marital_status_t1,
                health_t1,
                # sub-group interaction vars
                all_of(interactions_vars))

# check var structure is OK for glm
str(mi_subset2)

### recode vars for analysis
## outcomes
# change sf-12 vars to numeric
mi_subset2$sf12pcs_dv_t0 <- as.numeric(mi_subset2$sf12pcs_dv_t0)
mi_subset2$sf12mcs_dv_t0 <- as.numeric(mi_subset2$sf12mcs_dv_t0)
mi_subset2$sf12pcs_dv_t1 <- as.numeric(mi_subset2$sf12pcs_dv_t1)
mi_subset2$sf12mcs_dv_t1 <- as.numeric(mi_subset2$sf12mcs_dv_t1)

## make sure binary outcomes are in the right order
mi_subset2$srh_bin_t0 <- factor(mi_subset2$srh_bin_t0,
                             levels = c("excellent/very good", 
                                        "good/fair/poor"))
mi_subset2$srh_bin_t1 <- factor(mi_subset2$srh_bin_t1,
                             levels = c("excellent/very good", 
                                        "good/fair/poor"))

mi_subset2$ghq_case4_t0 <- factor(mi_subset2$ghq_case4_t0,
                               levels = c("0-3", "4 or more"))
mi_subset2$ghq_case4_t1 <- factor(mi_subset2$ghq_case4_t1,
                               levels = c("0-3", "4 or more"))

### change exposure to factor so it is ordered correctly for model
mi_subset2$exposure1 <- factor(mi_subset2$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))

### reorder covariates as required
## sex
mi_subset2$sex_dv_t0 <- factor(mi_subset2$sex_dv_t0,
levels = c("Female",
           "Male"))

## ethnicity
mi_subset2$non_white_t0 <- factor(mi_subset2$non_white_t0,
                               levels = c("White",
                                          "Non-white"))

## martial status
mi_subset2$marital_status_t0 <- factor(mi_subset2$marital_status_t0,
                               levels = c("married/civil partnership",
                                          "single",
                                          "divorced/separated/widowed"))

mi_subset2$marital_status_t1 <- factor(mi_subset2$marital_status_t1,
                                       levels = c("married/civil partnership",
                                                  "single",
                                                  "divorced/separated/widowed"))

## degree
mi_subset2$degree_bin_t0 <- factor(mi_subset2$degree_bin_t0,
                                   levels = c("No degree",
                                              "Degree or higher"))


## region
mi_subset2$gor_dv_t0 <- factor(mi_subset2$gor_dv_t0)

## industry
mi_subset2$sic2007_section_lab_t0 <- factor(mi_subset2$sic2007_section_lab_t0)

## occupation
mi_subset2$soc2000_major_group_title_t0 <- factor(mi_subset2$soc2000_major_group_title_t0)

## part-time/full-time
mi_subset2$jbft_dv_t0 <- factor(mi_subset2$jbft_dv_t0)#,
#                                levels = "FT employee", "PT employee")


## small firm
mi_subset2$small_firm_t0 <- factor(mi_subset2$small_firm_t0)#,
#                                levels = "over 50 employees under", "under 50 employees")


## employment contract
mi_subset2$emp_contract_t0 <- factor(mi_subset2$emp_contract_t0,
                                  levels = c("permanent",
                                             "fixed-term"))


## broken employment

mi_subset2$broken_emp_t0 <- factor(mi_subset2$broken_emp_t0,
                                     levels = c("Unbroken employment", "Broken employment"))

## second job
mi_subset2$j2has_dv_t0 <- factor(mi_subset2$j2has_dv_t0,
                                   levels = c("no", "yes"))

## relative poverty
mi_subset2$rel_pov_t0 <- factor(mi_subset2$rel_pov_t0,
                                 levels = c("not in relative poverty", "relative poverty"))

## long-term health condition
mi_subset2$health_t0 <- factor(mi_subset2$health_t0,
                               levels = c("no", "yes"))

mi_subset2$health_t1 <- factor(mi_subset2$health_t1,
                               levels = c("no", "yes"))

## reverse binary exposure var for matching 3:1
mi_subset2 <- mi_subset2 %>% 
  mutate(exp1_bin = ifelse(exposure1=="exposed (employed at t1)",
                           0,1)) # 1 = unexposed as in PS matching


# check var structure is OK for glm
str(mi_subset2)


mi_subset2_str <- mi_subset2 %>% 
  summary.default() %>% as.data.frame %>% 
  dplyr::group_by(Var1) %>%  
  tidyr::spread(key = Var2, value = Freq)

### take a random sample for testing code (if it's running slow)
#mi_subset2 <- sample_n(mi_subset2, 3000)

write_rds(mi_subset2, "./working_data/mi/mi_subset2.rds")

#### set method type depending on variable type --------------------------------
# specify method for each incomplete variable in data
# (numeric, binary, unordered, ordered)
# ("norm", "logreg", "polyreg", "polr")
mi_subset2_str$method <- ""
#mi_subset2_str$method[mi_subset2_str$Var1=="pidp"] <- ""
mi_subset2_str$method[mi_subset2_str$Var1=="sf12pcs_dv_t1"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="sf12mcs_dv_t1"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="srh_bin_t1"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="ghq_case4_t1"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="exposure1"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="sex_dv_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="age_dv_t0"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="non_white_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="marital_status_t0"] <- "polyreg"
mi_subset2_str$method[mi_subset2_str$Var1=="degree_bin_t0"] <- "logreg" 
mi_subset2_str$method[mi_subset2_str$Var1=="gor_dv_t0"] <- "polyreg"
mi_subset2_str$method[mi_subset2_str$Var1=="sic2007_section_lab_t0"] <- "polyreg"
mi_subset2_str$method[mi_subset2_str$Var1=="soc2000_major_group_title_t0"] <- "polyreg"
mi_subset2_str$method[mi_subset2_str$Var1=="jbft_dv_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="small_firm_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="emp_contract_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="broken_emp_t0"] <- "polr" 
mi_subset2_str$method[mi_subset2_str$Var1=="j2has_dv_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="rel_pov_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="health_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="sf12pcs_dv_t0"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="sf12mcs_dv_t0"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="srh_bin_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="ghq_case4_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="dep_child_bin_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="age_dv_t1"] <- "norm" # norm if including
mi_subset2_str$method[mi_subset2_str$Var1=="marital_status_t1"] <- "polyreg" #polyreg if including
mi_subset2_str$method[mi_subset2_str$Var1=="health_t1"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="exp1_bin"] <- ""
mi_subset2_str$method[mi_subset2_str$Var1=="sex_pcs"] <- "~I(sex_dv_t0*sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
mi_subset2_str$method[mi_subset2_str$Var1=="sex_mcs"] <- "~Isf12pcs_dv_t1-mean(sf12mcs_dv_t1, na.rm = TRUE)))"
mi_subset2_str$method[mi_subset2_str$Var1=="sex_srh"] <- "~I(sex_dv_t0*srh_bin2)"
mi_subset2_str$method[mi_subset2_str$Var1=="sex_ghq"] <- "~I(sex_dv_t0*ghq_bin)"
mi_subset2_str$method[mi_subset2_str$Var1=="age_pcs"] <- "~I(age_bin*sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
mi_subset2_str$method[mi_subset2_str$Var1=="age_mcs"] <- "~I(age_bin*sf12mcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
mi_subset2_str$method[mi_subset2_str$Var1=="age_srh"] <- "~I(age_bin*srh_bin2)"
mi_subset2_str$method[mi_subset2_str$Var1=="age_ghq"] <- "~I(age_bin*ghq_bin)"
mi_subset2_str$method[mi_subset2_str$Var1=="rel_pov_pcs"] <- "~I(rel_pov_t0*sf12pcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
mi_subset2_str$method[mi_subset2_str$Var1=="rel_pov_mcs"] <- "~I(rel_pov_t0*sf12mcs_dv_t1-mean(sf12pcs_dv_t1, na.rm = TRUE)))"
mi_subset2_str$method[mi_subset2_str$Var1=="rel_pov_srh"] <- "~I(rel_pov_t0*srh_bin2)"
mi_subset2_str$method[mi_subset2_str$Var1=="rel_pov_ghq"] <- "~I(rel_pov_t0*ghq_bin)"

## check mi_subset2_str order matches vars in mi_subset2
sum(as.vector(mi_subset2_str$Var1)!=as.vector(names(mi_subset2)))

## create vector to set new default methods methods for MI
myDefaultMethod <- as.vector(mi_subset2_str$method)

## define a custom predictorMatrix
# in matrix 1s are included in model; 0s are not
myPredictorMatrix <- make.predictorMatrix(mi_subset2)
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
myPredictorMatrix[c("sex_dv_t0","sf12pcs_dv_t1"),"sex_pcs"] <- 0
myPredictorMatrix[c("sex_dv_t0","sf12mcs_dv_t1"),"sex_pcs"] <- 0
myPredictorMatrix[c("sex_dv_t0","srh_bin_t1"),"sex_srh"] <- 0
myPredictorMatrix[c("sex_dv_t0","ghq_case4_t1"),"sex_ghq"] <- 0
myPredictorMatrix[c("age_dv_t0","sf12pcs_dv_t1"),"age_pcs"] <- 0
myPredictorMatrix[c("age_dv_t0","sf12mcs_dv_t1"),"age_pcs"] <- 0
myPredictorMatrix[c("age_dv_t0","srh_bin_t1"),"age_srh"] <- 0
myPredictorMatrix[c("age_dv_t0","ghq_case4_t1"),"age_ghq"] <- 0
myPredictorMatrix[c("rel_pov_t0","sf12pcs_dv_t1"),"rel_pov_pcs"] <- 0
myPredictorMatrix[c("rel_pov_t0","sf12mcs_dv_t1"),"rel_pov_pcs"] <- 0
myPredictorMatrix[c("rel_pov_t0","srh_bin_t1"),"rel_pov_srh"] <- 0
myPredictorMatrix[c("rel_pov_t0","ghq_case4_t1"),"rel_pov_ghq"] <- 0

myPredictorMatrix

#### imputation ----------------------------------------------------------------

set.seed(52267)
start_time <- Sys.time()
### for final MI model run 25 imputations and 10 iterations
## run fewer when testing code
imps2 <- mice(mi_subset2, m=25, maxit = 10,
#imps2 <- mice(mi_subset2, m=15, maxit = 10,
                            #defaultMethod=myDefaultMethod, 
              predictorMatrix=myPredictorMatrix,
             printFlag = FALSE)
end_time <- Sys.time()
end_time - start_time

summary(imps2)

### check NAs are removed
## convert data from mids format to standard df
temp2 <- complete(imps2, "long", include = FALSE)

## sum number of NAs by var
temp2_na <- temp2 %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x)))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(df="mi")

## compare with non-imputed data
mi_subset2_na <- mi_subset2  %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x)))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(df="og")

temp3 <- temp2_na %>% bind_rows(mi_subset2_na) %>% 
  pivot_wider(names_from = df, values_from = value)

write_csv(temp3, "./output/temp_output/temp3_micheck.csv")

## check NAs
sapply(complete(imps2,2), function(x) sum(is.na(x)))

#### save imputed data 
write_rds(imps2, "./working_data/mi/imputed_data.rds")



