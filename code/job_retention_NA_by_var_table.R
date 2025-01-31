################################################################################
######                   NAs by var supplementary table                   ######
################################################################################


library(tidyverse)

mi_subset2 <- readRDS("./working_data/mi/mi_subset2.rds")


mi_subset2_na <- mi_subset2  %>% 
  summarise(across(.cols=everything(),
                   .fns = ~sum(is.na(.x)))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(pc_na = (value/nrow(mi_subset2))*100)

write.csv(mi_subset2_na, "./output/cc/na_by_var.csv")
