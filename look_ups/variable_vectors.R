

#### vector for person identifier and weigting variables -----------------------
id_wt_vector <- c("pidp", "psu", "strata", "wt_name", "wt_value")

#### master vars vector (no time point suffix) ---------------------------------
#### vector baseline covariates ------------------------------------------------
master_var_vec <- c("sex_dv", 
                "age_dv",  
                "non_white", 
                "marital_status",
                "hiqual_dv", 
                "gor_dv",
                "sic2007_section_lab",
                "soc2000_major_group_title",
                "jbft_dv",
                "small_firm",
                "emp_contract",
                "broken_emp",
                "j2has_dv",
                "rel_pov",
                "health",
                "srh_bin",
                "ghq_case4",
                "sf12mcs_dv",
                "sf12pcs_dv",
                "jbstat", 
                "nunmpsp_dv")

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

#### vars including binary dummy vars for analysis -----------------------------
cov_vector3 <- c("sex_dv_t0" ,
                 "age_dv_t0" ,
                 "non_white_t0" ,
                 "marital_status_t0_married_civil_partnership" ,
                 "marital_status_t0_divorced_separated_widowed" ,
                 "marital_status_t0_single" ,
                 "hiqual_dv_t0_degree" ,
                 "hiqual_dv_t0_other_higher_degree" ,
                 "hiqual_dv_t0_a_level_etc" ,
                 "hiqual_dv_t0_gcse_etc" ,
                 "hiqual_dv_t0_other_qualification" ,
                 "hiqual_dv_t0_no_qualification" ,
                 "gor_dv_t0_east_midlands" ,
                 "gor_dv_t0_east_of_england" ,
                 "gor_dv_t0_london" ,
                 "gor_dv_t0_north_east" ,
                 "gor_dv_t0_north_west" ,
                 "gor_dv_t0_northern_ireland" ,
                 "gor_dv_t0_scotland" ,
                 "gor_dv_t0_south_east" ,
                 "gor_dv_t0_south_west" ,
                 "gor_dv_t0_wales" ,
                 "gor_dv_t0_west_midlands" ,
                 "gor_dv_t0_yorkshire_and_the_humber" ,
                 "sic2007_section_lab_t0_accommodation_and_food_service_activities" ,
                 "sic2007_section_lab_t0_administrative_and_support_service_activities" ,
                 "sic2007_section_lab_t0_construction" ,
                 "sic2007_section_lab_t0_education" ,
                 "sic2007_section_lab_t0_human_health_and_social_work_activities" ,
                 "sic2007_section_lab_t0_manufacturing" ,
                 "sic2007_section_lab_t0_other_industry" ,
                 "sic2007_section_lab_t0_professional_scientific_and_technical_activities" ,
                 "sic2007_section_lab_t0_public_administration_and_defence_compulsory_social_security" ,
                 "sic2007_section_lab_t0_transportation_and_storage" ,
                 "sic2007_section_lab_t0_wholesale_and_retail_trade_repair_of_motor_vehicles_and_motorcycles" ,
                 "soc2000_major_group_title_t0_administrative_and_secretarial_occupations" ,
                 "soc2000_major_group_title_t0_associate_professional_and_technical_occupations" ,
                 "soc2000_major_group_title_t0_elementary_occupations" ,
                 "soc2000_major_group_title_t0_managers_and_senior_officials" ,
                 "soc2000_major_group_title_t0_personal_service_occupations" ,
                 "soc2000_major_group_title_t0_process_plant_and_machine_operatives" ,
                 "soc2000_major_group_title_t0_sales_and_customer_service_occupations" ,
                 "soc2000_major_group_title_t0_science_and_technology_professionals" ,
                 "soc2000_major_group_title_t0_skilled_trades_occupations" ,
                 "jbft_dv_t0" ,
                 "small_firm_t0" ,
                 "emp_contract_t0" ,
                 "broken_emp_t0_broken_employment" ,
                 "broken_emp_t0_no_employment_spells" ,
                 "broken_emp_t0_unbroken_employment" ,
                 "j2has_dv_t0" ,
                 "rel_pov_t0" ,
                 "health_t0" ,
                 "srh_bin_t0" ,
                 "ghq_case4_t0" ,
                 "sf12mcs_dv_t0" ,
                 "sf12pcs_dv_t0")


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


