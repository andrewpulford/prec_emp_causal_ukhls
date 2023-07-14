################################################################################
# Function to create single digit variables for SIC 2007 and SOC 2000 
# classifications in UKHLS dataset




################################################################################
#####                          load look-up files                          #####
################################################################################

#### SIC 2007-------------------------------------------------------------------
sic2007 <- read.csv("C:/Users/0510028p/Documents/prec_emp_causal_ukhls/look_ups/publisheduksicsummaryofstructureworksheet.csv")

#### SOC 2000-------------------------------------------------------------------
soc2000 <- read.csv("C:/Users/0510028p/Documents/prec_emp_causal_ukhls/look_ups/soc2000_look-up.csv")


################################################################################
#####                          SIC 2007 data prep                         ######
################################################################################

### tidy things up
sic2007 <- sic2007 %>% clean_names() %>% 
  dplyr::select(description, section, level_headings) %>% 
  mutate(across(.cols = everything(), 
                .fns = ~str_to_sentence(.x))) %>% 
#    mutate(across(.cols = everything(),
#                .fns = ~gsub("[^ -~]", "", .x))) %>% 
    mutate(across(.cols = everything(),
                .fns = ~trimws(.x, which = "both"))) %>% 
  filter(level_headings%in%c("Section","Division"))

### temp df for sections
sic2007_section <- sic2007 %>% 
  filter(level_headings=="Section") %>% 
  dplyr::select(-level_headings) %>% 
  rename("sic2007_section_lab"="description")

### temp df for divisions
sic2007_division <- sic2007 %>% 
  filter(level_headings=="Division") %>% 
  dplyr::select(-level_headings) %>% 
  rename("jbsic07_cc"="description")

### join back together
sic2007 <- sic2007_division %>% 
  left_join(sic2007_section)

rm(sic2007_division, sic2007_section)

################################################################################
#####                          SOC 2000 data prep                         ######
################################################################################

### tidy things up
soc2000 <- soc2000 %>% clean_names() %>% 
  mutate(across(.cols = everything(), 
                .fns = ~str_to_sentence(.x))) %>% 
  mutate(across(.cols = everything(),
                .fns = ~trimws(.x, which = "both"))) %>% 
  rename("soc2000_major_group_title" = "major_group_title",
         "jbsoc00_cc" = "minor_group_title")

################################################################################
#####                          SIC 2007 function                          ######
################################################################################

sic2007f <- function(data = master_raw1){
  df <- data %>% 
    left_join(sic2007, by = "jbsic07_cc")  %>% 
    mutate(sic2007_section_lab = ifelse(jbsic07_cc %in% c("missing",                                                     
                                                  "inapplicable",                                                
                                                  "proxy",                                                       
                                                  "refusal",                                                     
                                                  "don't know",                                                  
                                                  "Only available for IEMB",                                     
                                                  "Not available for IEMB"),
                                   "missing", sic2007_section_lab)) %>% 
    mutate(sic2007_section_lab = replace_na(sic2007_section_lab, "missing"))
  
  
}

################################################################################
#####                          SOC 2000 function                          ######
################################################################################

soc2000f <- function(data = master_raw1){
  df <- data %>% 
    left_join(soc2000, by = "jbsoc00_cc")  %>% 
    mutate(soc2000_major_group_title = ifelse(jbsoc00_cc %in% c("missing",                                                     
                                                          "inapplicable",                                                
                                                          "proxy",                                                       
                                                          "refusal",                                                     
                                                          "don't know",                                                  
                                                          "Only available for IEMB",                                     
                                                          "Not available for IEMB"),
                                        "missing", soc2000_major_group_title)) %>% 
    mutate(soc2000_major_group_title = replace_na(soc2000_major_group_title, "missing"))
  
}

