################################################################################

# Persistent precarious employment and health - Understanding Society
# 1-02 - Create weight spine for causal analysis
# Andrew Pulford

# Data source:
# University of Essex, Institute for Social and Economic Research. (2021). 
# Understanding Society: Waves 1-10, 2009-2019 and Harmonised BHPS: Waves 1-18, 
# 1991-2009. [data collection]. 13th Edition. UK Data Service. SN: 6614, 
# http://doi.org/10.5255/UKDA-SN-6614-14

#### What this script does:
# (a) Creates weight spines for analytic samples a and b

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


### to be completed - refer back to descriptive weighting script
