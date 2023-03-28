################################################################################

# Precarious employment and health - Understanding Society
# 1-01 - Create master raw data and variable spine for causal analysis
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Creates a data dictionary of the variables to be included in the analytic dataset.  These were taken from the variable search section on the U-Soc website, and can be changed for other analyses if the column headings are retained.
# (b) Defines function "load_wave" to allow data frames from U-Soc waves to be created.  
# (c) Combines waves to create a master raw dataframe
# Data output: variable spine look-up, master raw dataframe

#### To be added:
# add analysis worksheet code for emp and inc
# "load_wave" calls for household responses

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
#####                         load and prepare data                        #####
################################################################################

var_path <- "C:/Users/0510028p/Documents/prec_emp_causal_ukhls/variables/"
data_path <- "C:/Users/0510028p/Documents/UKDA-6614-spss/spss/spss25/"


##----------------------------------------------------------------------------##
## Load data dictionary
##----------------------------------------------------------------------------##

vars <- read.csv(paste0(var_path,"analytic_vars.csv")) %>% 
  clean_names()

## save as RDF
write_rds(vars, "./variables/vars_master.rds")

##----------------------------------------------------------------------------##
## Load survey waves
##----------------------------------------------------------------------------##

#### define load_wave function
#     Function arguments:
#       wv_n - wave number
#       wv_l - wave letter
#       datafile - UNderstanding Society datafile

load_wave <- function(wv_n, wv_l, datafile) {
  temp_df <- read.spss(paste0(data_path,"ukhls_w",wv_n,"/",wv_l,"_indresp.sav"), 
                        to.data.frame=TRUE, use.value.labels=TRUE) %>% 
    as_tibble()

  ## select only vars in vars vector 
  
  # add wave prefix for all vars except pidp
  vars$variable <- if_else(vars$variable != "pidp", paste0(wv_l,"_",vars$variable), "pidp")
  # convert to vector
  var_selecter <- unlist(vars[,1])
  # create var names vector for wave
  names_use <- names(temp_df)[(names(temp_df) %in% var_selecter)]
  # select only vars to be used
  temp_df <- temp_df %>% dplyr::select(names_use) 
  
  # remove wave letter prefix from column names
  colnames(temp_df) <-sub ("^[^_]*_","",colnames(temp_df))
  
  # add vars for wave number and letter
  temp_df$wv_n <- wv_n
  temp_df$wv_l <- wv_l
  
  write_rds(temp_df,paste0("./raw_data/wv",wv_n,"_",datafile,"_raw.rds"))
}


#### Function calls
### Individual respondents
## wave 1 
load_wave(wv_n = 1, wv_l = "a", datafile = "indresp")

## wave 2  
load_wave(wv_n = 2, wv_l = "b", datafile = "indresp")

## wave 3
load_wave(wv_n = 3, wv_l = "c", datafile = "indresp")

## wave 4
load_wave(wv_n = 4, wv_l = "d", datafile = "indresp")

## wave 5
load_wave(wv_n = 5, wv_l = "e", datafile = "indresp")

## wave 6
load_wave(wv_n = 6, wv_l = "f", datafile = "indresp")

## wave 7
load_wave(wv_n = 7, wv_l = "g", datafile = "indresp")

## wave 8
load_wave(wv_n = 8, wv_l = "h", datafile = "indresp")

## wave 9
load_wave(wv_n = 9, wv_l = "i", datafile = "indresp")

## wave 10
load_wave(wv_n =10, wv_l = "j", datafile = "indresp")

##----------------------------------------------------------------------------##
## create single dataframe for waves 1-10 
##----------------------------------------------------------------------------##

## define list to add df's to
df_list <<- list()


load_raw_wave <- function(wv_n, wv_l, datafile){
  ## read file
  df <- readRDS(paste0("./raw_data/wv",wv_n,"_",datafile,"_raw.rds"))
  
  ## assign dataframes to global environment
  df_list[[wv_n]] <<- df
  
  #paste0("wv_",wv_n,"_raw")

  } # end of function --------------------

#### Individual respondents
# wave 1 function call ------------- 
load_raw_wave(wv_n = 1, wv_l = "a", datafile = "indresp")
# wave 2 function call ------------- 
load_raw_wave(wv_n = 2, wv_l = "b", datafile = "indresp")
# wave 3 function call -------------
load_raw_wave(wv_n = 3, wv_l = "c", datafile = "indresp")
# wave 4 function call -------------
load_raw_wave(wv_n = 4, wv_l = "d", datafile = "indresp")
# wave 5 function call -------------
load_raw_wave(wv_n = 5, wv_l = "e", datafile = "indresp")
# wave 6 function call -------------
load_raw_wave(wv_n = 6, wv_l = "f", datafile = "indresp")
# wave 7 function call -------------
load_raw_wave(wv_n = 7, wv_l = "g", datafile = "indresp")
# wave 8 function call -------------
load_raw_wave(wv_n = 8, wv_l = "h", datafile = "indresp")
# wave 9 function call -------------
load_raw_wave(wv_n = 9, wv_l = "i", datafile = "indresp")
# wave 10 function call -------------
load_raw_wave(wv_n = 10, wv_l = "j", datafile = "indresp")


## bind listed dfs into single master raw df
master_raw1 <- bind_rows(df_list)#, .id = "column_label")



## write as RDS
write_rds(master_raw1, "./raw_data/master_raw1.rds")

# remove list
rm(df_list)

### Households
## to be added?





