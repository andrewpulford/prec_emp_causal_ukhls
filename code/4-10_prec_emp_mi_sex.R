################################################################################

# Precarious employment and health - Understanding Society
# 4-04 - IPTW multiple imputed analytic sample - sub-group analysis by sex
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

################################################################################
#####                         load and prepare data                        #####
################################################################################


#### read in variable vectors --------------------------------------------------
source("./look_ups/variable_vectors.r")

##### load original data without imoutation ------------------------------------
#mi_subset2 <-  readRDS("./working_data/mi/mi_subset2.rds")

#### load imputed data --------------------------------------------------------
imputed_data <- readRDS("./working_data/mi/imputed_data.rds")


#### prepare data -------------------------------------------------------------- 

## create male and female dfs
df_long <- complete(imputed_data, "long",include = T)

sapply(df_long, function(x) sum(is.na(x)))


df_long_f <- df_long[which(df_long$sex_dv_t0 == 'Female'),]
df_long_m <- df_long[which(df_long$sex_dv_t0 == 'Male'),]

imp_f <- as.mids(df_long_f)
imp_m <- as.mids(df_long_m)

sapply(complete(imp_f,"long"), function(x) sum(is.na(x)))
sapply(complete(imp_m,"long"), function(x) sum(is.na(x)))


################################################################################
#####               inverse probability of treatment weighting             #####
################################################################################

#### function ------------------------------------------------------------------

iptw_func_sex <- function(data){
                   weightthem(exp1_bin ~
 #                           sex_dv_t0 +
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
                            rel_pov_t0 +
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

weightit_f <- iptw_func_sex(data = imp_f)
summary(weightit_f)

#with(weightit_df, summary(as.data.frame(mget(ls()))))

#imputed_data$weights_ps <- weightit_df$weights

#with(imputed_data, sum(imputed_data$weights_ps))

test <- bal.tab(weightit_f, un = TRUE, 
                binary = "std", continuous = "std")
test2 <- test$Balance.Across.Imputations

## probably don't need these....
#bal.plot(weightit_f, which.imp = 1, 
#         var.name = "sex_dv_t0", 
#         which = "both")
bal.plot(weightit_f, which.imp = 1, 
         var.name = "age_dv_t0", 
         which = "both")
bal.plot(weightit_f, which.imp = 1, 
         var.name = "non_white_t0", 
         which = "both")
bal.plot(weightit_f, which.imp = 1, 
         var.name = "rel_pov_t0", 
         which = "both")

## create love plot to visualise balance between unmatched and matched data across MIs
tiff("./output/mi/sub_groups/sex/f_love_plot.tiff", width = 960)
love.plot(weightit_f, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
dev.off()
# var.names() - to clean up names
# https://www.rdocumentation.org/packages/cobalt/versions/3.4.1/topics/var.names

#### males -------------------------------------------------------------------

weightit_m <- iptw_func_sex(data = imp_m)
summary(weightit_m)

#with(weightit_df, summary(as.data.frame(mget(ls()))))

#imputed_data$weights_ps <- weightit_df$weights

#with(imputed_data, sum(imputed_data$weights_ps))

test <- bal.tab(weightit_m, un = TRUE, 
                binary = "std", continuous = "std")
test2 <- test$Balance.Across.Imputations

## probably don't need these....
#bal.plot(weightit_f, which.imp = 1, 
#         var.name = "sex_dv_t0", 
#         which = "both")
bal.plot(weightit_m, which.imp = 1, 
         var.name = "age_dv_t0", 
         which = "both")
bal.plot(weightit_m, which.imp = 1, 
         var.name = "non_white_t0", 
         which = "both")
bal.plot(weightit_m, which.imp = 1, 
         var.name = "rel_pov_t0", 
         which = "both")

## create love plot to visualise balance between unmatched and matched data across MIs
tiff("./output/mi/sub_groups/sex/m_love_plot.tiff", width = 960)
love.plot(weightit_m, thresholds = 0.1, stats = "m",
          drop.distance = TRUE, binary = "std", continuous = "std",
          sample.names = c("Unweighted","Inverse probability of treatment weighted")) 
dev.off()
# var.names() - to clean up names
# https://www.rdocumentation.org/packages/cobalt/versions/3.4.1/topics/var.names

################################################################################
#####                               Descriptives                           #####
################################################################################

### leave for now
# use createtableone and leave SMDs for next script
#table_one_ps <- with(weightit_df, CreateTableOne(
#  vars = cov_vector,
#  data=as.data.frame(mget(ls())), 
#  factorVars=catVars_short_vec,
#  strata = "exposure1", test =TRUE))


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
#                                    sex_dv_t0 +
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
                                    rel_pov_t0 +
                                    health_t0 +
                                    health_t1 +
                                    # interaction terms
#                                    sex_dv_t0*age_dv_t0 +
#                                    sex_dv_t0*rel_pov_t0 +
                                    age_dv_t0*rel_pov_t0 +
                                    (1|pidp)))
}

### sf-12 mcs
iptw_dr_mcs <- function(data){
  with(data,
       glmmTMB(sf12mcs_dv_t1 ~
                 exposure1 +
                 sf12mcs_dv_t0 +
                 #                                    sex_dv_t0 +
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
                 rel_pov_t0 +
                 health_t0 +
                 health_t1 +
                 # interaction terms
                 #                                    sex_dv_t0*age_dv_t0 +
                 #                                    sex_dv_t0*rel_pov_t0 +
                 age_dv_t0*rel_pov_t0 +
                 (1|pidp)))
}

### self-rated health
iptw_dr_srh <- function(data){
  with(data,
       glmmTMB(srh_bin_t1 ~
                 exposure1 +
                 srh_bin_t0 +
                 #                                    sex_dv_t0 +
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
                 rel_pov_t0 +
                 health_t0 +
                 health_t1 +
                 # interaction terms
                 #                                    sex_dv_t0*age_dv_t0 +
                 #                                    sex_dv_t0*rel_pov_t0 +
                 age_dv_t0*rel_pov_t0 +
                 (1|pidp)))
}

### ghq-12
iptw_dr_ghq <- function(data){
  with(data,
       glmmTMB(ghq_case4_t1 ~
                 exposure1 +
                 ghq_case4_t0 +
                 #                                    sex_dv_t0 +
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
                 rel_pov_t0 +
                 health_t0 +
                 health_t1 +
                 # interaction terms
                 #                                    sex_dv_t0*age_dv_t0 +
                 #                                    sex_dv_t0*rel_pov_t0 +
                 age_dv_t0*rel_pov_t0 +
                 (1|pidp)))
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

#here <<<<<<<<<<<<<<<<<<< add ocde for males


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
  dplyr::select(c(outcome,	est_type,	estimate,	std.error,	p.value,	lci,	uci))

write.csv(combined_df, "./output/mi/weighted_outcomes/mi_dr_iptw_df.csv")


##################


sapply(complete(imputed_data,"long"), function(x) sum(is.na(x)))
