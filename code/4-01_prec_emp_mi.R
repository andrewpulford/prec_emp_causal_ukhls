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

################################################################################
#####                         load and prepare data                        #####
################################################################################


####read in variable vectors ---------------------------------------------------
source("./look_ups/variable_vectors.r")

## extra vars for creating exposure vars
extra_vars <- c("jbstat_t1", "nunmpsp_dv_t1")



####load eligible cases --------------------------------------------------------
pair_eligible <- readRDS("./working_data/pair_eligible.rds") %>% 
    dplyr::select(pidp, exposure1, all_of(c(cov_vector, cov_vector2, 
                           outcome_vector2)))



# check number of missing cases by vartiable
sapply(pair_eligible, function(x) sum(is.na(x)))

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
            "nunmpsp_dv_t0")

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




################################################################################
#####                          Multiple imputation                         #####
################################################################################

#### MI for one variable -------------------------------------------------------

### create sub-set of variables and cases
mi_subset1 <- pair_eligible %>% dplyr::select(sf12pcs_dv_t1, sex_dv_t0, age_dv_t0, exposure1)

# check var structure is OK for glm
str(mi_subset1)

# change outcome var to numeric
mi_subset1$sf12pcs_dv_t1 <- as.numeric(mi_subset1$sf12pcs_dv_t1)
# change exposure to factor so it is ordered correctly for model
mi_subset1$exposure1 <- factor(mi_subset1$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))

# take a random sample for testing code (if it's running slow)
#mi_subset1 <- sample_n(mi_subset1, 1000)

## set seed
set.seed(7341)

## create imputations
imps <- mice(mi_subset1, m=10, maxit=1)
summary(imps)

### Analysing imputed datasets

## log reg on SF-12 PCS with one covariate
fit <- with(imps, exp=glm(sf12pcs_dv_t1~sex_dv_t0, family="gaussian"))
summary(pool(fit), conf.int = TRUE)

#### MI for multiple variables --------------------------------------------------------

### create sub-set of variables and cases
mi_subset2 <- pair_eligible %>% 
  dplyr::select(
                pidp,
                sf12pcs_dv_t1, 
                exposure1, 
                all_of(cov_vector),
                age_dv_t1,
                marital_status_t1,
                gor_dv_t1,
                health_t1,
                srh_bin_t1,
                ghq_case4_t1,
                sf12mcs_dv_t1
                )

# check var structure is OK for glm
str(mi_subset2)

# change sf-12 vars to numeric
mi_subset2$sf12pcs_dv_t0 <- as.numeric(mi_subset2$sf12pcs_dv_t0)
mi_subset2$sf12mcs_dv_t0 <- as.numeric(mi_subset2$sf12mcs_dv_t0)
mi_subset2$sf12pcs_dv_t1 <- as.numeric(mi_subset2$sf12pcs_dv_t1)
mi_subset2$sf12mcs_dv_t1 <- as.numeric(mi_subset2$sf12mcs_dv_t1)
# change exposure to factor so it is ordered correctly for model
mi_subset2$exposure1 <- factor(mi_subset2$exposure1,
                               levels = c("unexposed",
                                          "exposed (employed at t1)"))

# take a random sample for testing code (if it's running slow)
mi_subset2 <- sample_n(mi_subset2, 3000)

mi_subset2_str <- mi_subset2 %>% 
  summary.default() %>% as.data.frame %>% 
  dplyr::group_by(Var1) %>%  
  tidyr::spread(key = Var2, value = Freq)

## set method type depending on variable type
mi_subset2_str$method <- ""
#mi_subset2_str$method[mi_subset2_str$Var1=="pidp"] <- ""
mi_subset2_str$method[mi_subset2_str$Var1=="sf12pcs_dv_t1"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="exposure1"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="sex_dv_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="age_dv_t0"] <- "norm"
mi_subset2_str$method[mi_subset2_str$Var1=="non_white_t0"] <- "logreg"
mi_subset2_str$method[mi_subset2_str$Var1=="marital_status_t0"] <- "polyreg"
mi_subset2_str$method[mi_subset2_str$Var1=="hiqual_dv_t0"] <- "polr"
mi_subset2_str$method[mi_subset2_str$Var1=="gor_dv_t0"] <- "polyreg"


## set new default methods
# specify method for each incomplete variable in data
# (numeric, binary, unordered, ordered)
myDefaultMethod <- c("norm", "logreg", "polyreg", "polr")

## define a custom predictorMatrix
# in matrix 1s are included in model; 0s are not
myPredictorMatrix <- make.predictorMatrix(mi_subset2)
myPredictorMatrix[,"pidp"] <- 0
myPredictorMatrix["pidp",] <- 0
myPredictorMatrix

## We can now proceed to imputation, being careful to pass myDefaultMethod and 
## myPredictorMatrix to the mice function:

set.seed(52267)
imps2 <- mice(mi_subset2, m=5, maxit = 10,
              defaultMethod=myDefaultMethod, predictorMatrix=myPredictorMatrix,
             printFlag = FALSE)
summary(imps2)

## We can then fit our substantive model to the imputations and pool the results:
fit <- with(data = imps2, exp = glmmTMB(sf12pcs_dv_t1 ~
                                     exposure1 +
                                     sex_dv_t0 +
                                     age_dv_t0 +
                                     non_white_t0 +
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
                                     age_dv_t1 +
                                     marital_status_t1 +
                                     gor_dv_t1 +
                                     health_t1 +
                                     srh_bin_t1 +
                                     ghq_case4_t1 +
                                     sf12mcs_dv_t1 +
                                     # interaction terms
                                     sex_dv_t0*age_dv_t0 +
                                     sex_dv_t0*rel_pov_t0 +
                                     age_dv_t0*rel_pov_t0 +
                                     (1|pidp), 
                                   family="gaussian"))
pooled <- pool(fit)
mi_mod2 <- summary(pooled, conf.int = TRUE)

# add pidp back into matrix for mlm - set pidp to zero?
# number of MIs and iterations?
