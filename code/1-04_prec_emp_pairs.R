################################################################################

# Precarious employment and health - Understanding Society
# 1-04 - Create paired and three-wave clean master raw dataframes for causal analysis
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) creates potential person-trials dataframes.
# (b) .  
# (c) .
# Data output: paired and three-wave cleaned master raw dataframe(s)

#### To be added:

################################################################################


## remove any existing objects from global environment
rm(list=ls()) 


################################################################################
#####                            install packages                          #####
################################################################################

library(tidyverse) # all kinds of stuff 
library(janitor) # cleaning up


################################################################################
#####                               load data                              #####
################################################################################

clean_raw1 <- readRDS("./working_data/master_raw1_clean.rds")

#clean_raw1 %>% filter(wv_n == 1 | wv_n ==2)

################################################################################
#####                           define functions                           #####
################################################################################

#### ---------------------------------------------------------------------------
#### Paired waves
#### ---------------------------------------------------------------------------


#### creates raw df for paired cases -------------------------------------------
# includes participants with wave response at either t0 or t1
pair_func <- function(t0, t1){
  df <- clean_raw1 %>% filter(wv_n == t0 | wv_n == t1) %>% 
    arrange(pidp,wv_n) %>% 
    group_by(pidp) %>%
    mutate(t0_flag = ifelse(wv_n == t0,1,0),
           t1_flag = ifelse(wv_n == t1,1,0)) %>% 
    mutate(n_rows = n()) %>% # create var indicating number of rows per pidp
    ungroup() 
}

#### creates wide version raw paired cases -------------------------------------
# note currently uses left join so only valid t0 cases are included
# change to full_join if t0 missing cases to be included
pair_wide_func <- function(data){
  df1 <- data %>% filter(t0_flag==1) %>% 
    rename_at(vars(-c(pidp, t0_flag, t1_flag)), ~ paste0(.,"_t0"))
  
  df2 <- data %>% filter(t1_flag==1) %>% 
    rename_at(vars(-c(pidp, t0_flag, t1_flag)), ~ paste0(.,"_t1"))
  
  df <- df1 %>% left_join(df2) %>%
    mutate(t0_flag = ifelse(is.na(t0_flag),0,t0_flag),
           t1_flag = ifelse(is.na(t1_flag),0,t1_flag)) %>% 
    mutate(t_flag = paste0(t0_flag,"-",t1_flag))
  
}

#### creates raw df for paired complete case analysis --------------------------
# includes participants with wave response at both t0 and t1
pair_cc_func <- function(data, pair){
  df <- data %>% 
    filter(n_rows==2) %>% 
    dplyr::select(-n_rows) %>% 
    mutate(wave_pair = pair) %>% 
    group_by(pidp) %>% 
    mutate(t = row_number(),
           t = ifelse(t==1,"t0",
                      ifelse(t==2,"t1",
                             "check"))) %>% 
    ungroup()
}

#### creates wide version of complete case data --------------------------------
pair_cc_wide_func <- function(data){
  df1 <- data %>% filter(t=="t0") %>% 
    rename_at(vars(-pidp), ~ paste0(.,"_t0"))
  
  df2 <- data %>% filter(t=="t1") %>% 
    rename_at(vars(-pidp), ~ paste0(.,"_t1"))
  
  df <- df1 %>% left_join(df2)
  
  }

#### ---------------------------------------------------------------------------
#### Three waves
#### ---------------------------------------------------------------------------

#### creates three wave raw dataframes -----------------------------------------
# includes participants with wave response at either t0 or t1
three_func <- function(t0, t1, t2){
  df <- clean_raw1 %>% filter(wv_n == t0 | wv_n == t1 | wv_n == t2) %>% 
    arrange(pidp,wv_n) %>% 
    group_by(pidp) %>%
    mutate(t0_flag = ifelse(wv_n == t0,1,0),
           t1_flag = ifelse(wv_n == t1,1,0),
           t2_flag = ifelse(wv_n == t2,1,0)) %>% 
    mutate(n_rows = n()) %>% # create var indicating number of rows per pidp
    ungroup() 
}



#### creates wide three wave raw dataframes ------------------------------------
# note currently uses left join so only valid t0 cases are included
# change to full_join if t0 missing cases to be included
three_wide_func <- function(data){
  df1 <- data %>% filter(t0_flag==1) %>% 
    rename_at(vars(-c(pidp, t0_flag, t1_flag, t2_flag)), ~ paste0(.,"_t0"))
  
  df2 <- data %>% filter(t1_flag==1) %>% 
    rename_at(vars(-c(pidp, t0_flag, t1_flag, t2_flag)), ~ paste0(.,"_t1"))
  
  df3 <- data %>% filter(t2_flag==1) %>% 
    rename_at(vars(-c(pidp, t0_flag, t1_flag, t2_flag)), ~ paste0(.,"_t2"))
  
  df <- df1 %>% left_join(df2) %>%
    mutate(t0_flag = ifelse(is.na(t0_flag),0,t0_flag),
           t1_flag = ifelse(is.na(t1_flag),0,t1_flag)) %>% 
    mutate(t_flag = paste0(t0_flag,"-",t1_flag))
  
}


#### creates three wave complete case dataframes -------------------------------
three_cc_func <- function(data, trio){
  df <- data %>% 
    filter(n_rows==3) %>% 
    dplyr::select(-n_rows) %>% 
    mutate(wave_trio = trio) %>% 
    group_by(pidp) %>% 
    mutate(t = row_number(),
           t = ifelse(t==1,"t0",
                      ifelse(t==2,"t1",
                             ifelse(t==3,"t2",
                             "check")))) %>% 
    ungroup()
}


#### creates wide three wave complete case dataframes --------------------------
three_cc_wide_func <- function(data){
  df1 <- data %>% filter(t=="t0") %>% 
    rename_at(vars(-pidp), ~ paste0(.,"_t0"))
  
  df2 <- data %>% filter(t=="t1") %>% 
    rename_at(vars(-pidp), ~ paste0(.,"_t1"))
  
  df3 <- data %>% filter(t=="t2") %>% 
    rename_at(vars(-pidp), ~ paste0(.,"_t2"))
  
  df <- df1 %>% left_join(df2, df3, by = "pidp")
  
}


################################################################################
#####                            function calls                            #####
################################################################################

#### ---------------------------------------------------------------------------
#### Paired waves
#### ---------------------------------------------------------------------------

#### waves 1 and 2 -------------------------------------------------------------
#pair12_raw <- pair_func(t0 = 1,t1 = 2)
#pair12_raw_wide <- pair_wide_func(data=pair12_raw)
#pair12_cc <- pair_cc_func(data = pair12_raw, pair = "pair12")
#pair12_cc_wide <- pair_cc_wide_func(data=pair12_cc)

#### waves 2 and 3 -------------------------------------------------------------
pair23_raw <- pair_func(t0 = 2,t1 = 3)
pair23_raw_wide <- pair_wide_func(data=pair23_raw)
pair23_cc <- pair_cc_func(data = pair23_raw, pair = "pair23")
pair23_cc_wide <- pair_cc_wide_func(data=pair23_cc)

#### waves 3 and 4 -------------------------------------------------------------
pair34_raw <- pair_func(t0 = 3,t1 = 4)
pair34_raw_wide <- pair_wide_func(data=pair34_raw)
pair34_cc <- pair_cc_func(data = pair34_raw, pair = "pair34")
pair34_cc_wide <- pair_cc_wide_func(data=pair34_cc)

#### waves 4 and 5 -------------------------------------------------------------
pair45_raw <- pair_func(t0 = 4,t1 = 5)
pair45_raw_wide <- pair_wide_func(data=pair45_raw)
pair45_cc <- pair_cc_func(data = pair45_raw, pair = "pair45")
pair45_cc_wide <- pair_cc_wide_func(data=pair45_cc)

#### waves 5 and 6 -------------------------------------------------------------
pair56_raw <- pair_func(t0 = 5,t1 = 6)
pair56_raw_wide <- pair_wide_func(data=pair56_raw)
pair56_cc <- pair_cc_func(data = pair56_raw, pair = "pair56")
pair56_cc_wide <- pair_cc_wide_func(data=pair56_cc)

#### waves 6 and 7 -------------------------------------------------------------
pair67_raw <- pair_func(t0 = 6,t1 = 7)
pair67_raw_wide <- pair_wide_func(data=pair67_raw)
pair67_cc <- pair_cc_func(data = pair67_raw, pair = "pair67")
pair67_cc_wide <- pair_cc_wide_func(data=pair67_cc)

#### waves 7 and 8 -------------------------------------------------------------
pair78_raw <- pair_func(t0 = 7,t1 = 8)
pair78_raw_wide <- pair_wide_func(data=pair78_raw)
pair78_cc <- pair_cc_func(data = pair78_raw, pair = "pair78")
pair78_cc_wide <- pair_cc_wide_func(data=pair78_cc)

#### waves 8 and 9 -------------------------------------------------------------
pair89_raw <- pair_func(t0 = 8,t1 = 9)
pair89_raw_wide <- pair_wide_func(data=pair89_raw)
pair89_cc <- pair_cc_func(data = pair89_raw, pair = "pair89")
pair89_cc_wide <- pair_cc_wide_func(data=pair89_cc)

#### waves 9 and 10 ------------------------------------------------------------
pair910_raw <- pair_func(t0 = 9,t1 = 10)
pair910_raw_wide <- pair_wide_func(data=pair910_raw)
pair910_cc <- pair_cc_func(data = pair910_raw, pair = "pair910")
pair910_cc_wide <- pair_cc_wide_func(data=pair910_cc)

#### ---------------------------------------------------------------------------
#### Three waves
#### ---------------------------------------------------------------------------
#### waves 1, 2 and 3 ----------------------------------------------------------
#three123_raw <- three_func(t0=1, t1=2, t2=3)
#three123_wide_raw <- three_wide_func(data=three123_raw)
#three123_cc <- three_cc_func(data=three123_raw, trio = "trio123")
#three123_cc_wide <- three_cc_wide_func(three123_cc)

#### waves 2, 3 and 4 ----------------------------------------------------------
three234_raw <- three_func(t0=2, t1=3, t2=4)
three234_wide_raw <- three_wide_func(data=three234_raw)
three234_cc <- three_cc_func(data=three234_raw, trio = "trio234")
three234_cc_wide <- three_cc_wide_func(three234_cc)

#### waves 3, 4 and 5 ----------------------------------------------------------
three345_raw <- three_func(t0=3, t1=4, t2=5)
three345_wide_raw <- three_wide_func(data=three345_raw)
three345_cc <- three_cc_func(data=three345_raw, trio = "trio345")
three345_cc_wide <- three_cc_wide_func(three345_cc)

#### waves 4, 5 and 6 ----------------------------------------------------------
three456_raw <- three_func(t0=4, t1=5, t2=6)
three456_wide_raw <- three_wide_func(data=three456_raw)
three456_cc <- three_cc_func(data=three456_raw, trio = "trio456")
three456_cc_wide <- three_cc_wide_func(three456_cc)

#### waves 5, 6 and 7 ----------------------------------------------------------
three567_raw <- three_func(t0=5, t1=6, t2=7)
three567_wide_raw <- three_wide_func(data=three567_raw)
three567_cc <- three_cc_func(data=three567_raw, trio = "trio567")
three567_cc_wide <- three_cc_wide_func(three567_cc)

#### waves 6, 7 and 8 ----------------------------------------------------------
three678_raw <- three_func(t0=6, t1=7, t2=8)
three678_wide_raw <- three_wide_func(data=three678_raw)
three678_cc <- three_cc_func(data=three678_raw, trio = "trio678")
three678_cc_wide <- three_cc_wide_func(three678_cc)

#### waves 7, 8 and 9 ----------------------------------------------------------
three789_raw <- three_func(t0=7, t1=8, t2=9)
three789_wide_raw <- three_wide_func(data=three789_raw)
three789_cc <- three_cc_func(data=three789_raw, trio = "trio789")
three789_cc_wide <- three_cc_wide_func(three789_cc)

#### waves 8, 9 and 10 ---------------------------------------------------------
three8910_raw <- three_func(t0=8, t1=9, t2=10)
three8910_wide_raw <- three_wide_func(data=three8910_raw)
three8910_cc <- three_cc_func(data=three8910_raw, trio = "trio8910")
three8910_cc_wide <- three_cc_wide_func(three8910_cc)



################################################################################
#####                 Create paired complete cases dataset                 #####
################################################################################

### bind into single df
pair_cc_raw <- pair23_cc_wide %>% 
  bind_rows(pair34_cc_wide,
            pair45_cc_wide,
            pair56_cc_wide,
            pair67_cc_wide,
            pair78_cc_wide,
            pair89_cc_wide,
            pair910_cc_wide)

# save df
write_rds(pair_cc_raw, "./working_data/pair_cc_raw.rds")

### create list of df's
#list()
# leave unless needed
################################################################################
#####                  Create trio complete cases dataset                  #####
################################################################################

three_cc_raw <- three234_cc_wide %>% 
  bind_rows(three345_cc_wide,
            three456_cc_wide,
            three567_cc_wide,
            three678_cc_wide,
            three789_cc_wide,
            three8910_cc_wide)
