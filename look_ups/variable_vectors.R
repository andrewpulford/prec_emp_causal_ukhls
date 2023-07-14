

#### vector for person identifier and weigting variables -----------------------
id_wt_vector <- c("pidp", "psu", "strata", "wt_name", "wt_value")

#### vector baseline covariates ------------------------------------------------
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

#### vector vor covariates measured at t1 --------------------------------------
cov_vector2 <- c("age_dv_t1",  
                 "marital_status_t1",
                 "gor_dv_t1",
                 "health_t1",
                 "srh_bin_t1",
                 "ghq_case4_t1",
                 "sf12mcs_dv_t1",
                 "sf12pcs_dv_t1")

#### vector for outcome analysis (inc exposure variables) ----------------------
outcome_vector <- c("sf12pcs_dv_t1",
                     "sf12mcs_dv_t1",
                     "srh_bin_t1",
                     "ghq_case4_t1",
                    "exposure1",
                    "exposure2")

#### vector for outcome variables only------------------------------------------
outcome_vector2 <- c("sf12pcs_dv_t1",
                     "sf12mcs_dv_t1",
                     "srh_bin_t1",
                     "ghq_case4_t1")

#######################xxxxxxxxxxxxxxxxxx


