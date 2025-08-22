#### ---------------------------------------------------------------------- ####
###                                   DAG                                    ###
#### ---------------------------------------------------------------------- ####

## remove any existing objects from global environment
rm(list=ls()) 

## disable scientific notation printing
options(scipen=999)

## load packages
library(DiagrammeR)
library(tidyverse)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)
library(xml2)



################################################################################
#####                              create DAG                              #####
################################################################################


#### Mental health outcomes ----------------------------------------------------

DAGJRS <-DiagrammeR::grViz("
                  digraph{
                  graph[ranksep=0.2]
                  
                  node[shape=plaintext,fontname=Arial]
                  
                    Inv[label ='Time invariant confounders',shape=oval]
                    Var0[label='Time-varying \nconfounders T0',shape=box]
                    Var1[label='Time-varying \nconfounders T1',shape=box]
                    Age0[label='Age T0', shape=oval]
                    Age1[label='Age T1', shape=oval]
                    Emp0[label='Employment status \nat T0']
                    JRS[label=<<B>Job retention scheme</B>>,fontcolor=darkgreen]
                    MH0[label='Mental Health T0']
                    MH1[label=<<B>Mental Health T1</B>>,fontcolor=darkorange2]

                  edge[minlen=2]
                    Inv->Var0->MH1
                    Inv->Var1
                    Inv->MH0->MH1
                    Inv->MH1
                    Inv->Emp0
                    Inv->JRS
                    Var0->Var1
                    Var0->MH1
                    Var0->Emp0
                    Var1->JRS
                    MH0->Var0
                    Var1->MH1
                    Emp0->MH1
                    Emp0->JRS
                    MH0->Emp0
                    JRS->MH1
                    Age0->Emp0
                    Age0->MH0
                    Age0->Var0
                    Age1->JRS
                    Age1->MH1
                    Age1->Var1
                    Emp0->Var1
                    
                    
                    
                  {rank=min; Inv}
                  {rank=same; Var0; Var1; Age0; Age1}
                  {rank=same; Emp0; JRS}
                  {rank=max; MH0; MH1}
                  }
                  ")

DAGJRS

DAGJRS %>% 
  export_svg() %>% 
  charToRaw() %>% 
  rsvg_png("./output/job_retention_tte_dag_mh.png")

#### Alt version ---------------------------------------------------------------

DAGJRS2 <-DiagrammeR::grViz("
                  digraph{
                  graph[ranksep=0.2]
                  
                  node[shape=plaintext,fontname=Arial]
                  
                    Inv[label ='Time invariant confounders',shape=oval]
                    Var0[label='Time-varying \nconfounders T0',shape=box]
                    Var1[label='Time-varying \nconfounders T1',shape=box]
                    Age0[label='Age T0', shape=oval]
                    Age1[label='Age T1', shape=oval]
                    Emp0[label='Employment status \nat T0']
                    JRS[label=<<B>Job retention scheme</B>>,fontcolor=darkgreen]
                    MH0[label='Mental Health T0']
                    MH1[label=<<B>Mental Health T1</B>>,fontcolor=darkorange2]

                  edge[minlen=2]
                    Inv->Var0->MH1
                    Inv->Var1
                    Inv->MH0->MH1
                    Inv->MH1
                    Inv->Emp0
                    Var0->Var1
                    Var0->MH1
                    Var0->Emp0
                    MH0->Var0
                    Var1->MH1
                    Emp0->MH1
                    Emp0->JRS
                    MH0->Emp0
                    JRS->MH1
                    Age0->Emp0
                    Age0->MH0
                    Age0->Var0
                    Age1->MH1
                    Age1->Var1
                    Emp0->Var1
                    
                    
                    
                  {rank=min; Inv}
                  {rank=same; Var0; Var1; Age0; Age1}
                  {rank=same; Emp0; JRS}
                  {rank=max; MH0; MH1}
                  }
                  ") 

DAGJRS2

DAGJRS2 %>%
  export_svg() %>% 
  charToRaw() %>% 
  rsvg_png("./output/job_retention_tte_dag_mh_ALT.png")

#### DAG key -------------------------------------------------------------------
DAGJRS_key <- DiagrammeR::grViz("
                  digraph{
                  graph[ranksep=0]
                  
                  node[shape=plaintext,fontname=Arial]
                  
                  Key1[label='Time invariant confounders: baseline outcome measure, gender, ethnicity,\neducational attainment, industry classification,\noccupational classification, size of workplace, full or part-time working,\nusual working hours, disability or long-term condition, employment contract,\nemployment continuity since previous wave,\nmultiple employment and relative poverty (<60% of median income)',shape=box]
                  Key2[label='Time-varying covariates: age, marital status, region and dependent children',shape=box]
                  
                  edge[minlen=2]
                  Key1->Key2 [style=invis]
                  
                  {rank=max; Key2}
                  }
                  ")

DAGJRS_key

DAGJRS_key %>%
  export_svg() %>% 
  charToRaw() %>% 
  rsvg_png("./output/job_retention_tte_dag_mh_key.png")


