#### ---------------------------------------------------------------------- ####
###                              STROBE diagram                              ###
#### ---------------------------------------------------------------------- ####

## remove any existing objects from global environment
rm(list=ls()) 

## disable scientific notation printing
options(scipen=999)

## load packages
library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)


################################################################################
#### create flowchart
strobe <- grViz("
digraph {
  graph [ranksep = 0.2]

  node [shape = box, width = 3, fontname = arial, fontsize = 12]
    F [label = '29,456 lost to follow-up:\\n
    3,584 self-employed or undertaking\\nunpaid work for family business at t1\\n
    134	government training scheme or\\napprenticeship at t1\\n
    1,316	full-time student at t1\\n
    1,936	retired at t1\\n
    1,328	on maternity leave at t1\\n
    974	caring for family or home at t1\\n
    518	long-term sick or disabled at t1\\n
    279	other unspecified employment\\nstatus at t1\\n
    19,387	health outcomes not measured\\nat t0 or t1
    ']
    A [label = '300,409 potential person-trials screened for eligibility']
    B [label = '81,741 potential person-trials not eligible for trial:\\n
    71,960 not working-age\\n
    91,522 not employee at t0']
    C [label = '136,927 eligible potential person-trials']
    D [label = '134,243 allocated to treatment group']
    E [label = '2,684 allocated to control group']
    
    G [label = '641 lost to follow-up:\\n
    5	retired at t1\\n
    636	health outcomes not measured at t0 or t1']
    H[label = '104,787 potential-person trial outcomes assessed']
    I[label = '2,043 potential-person trial outcomes assessed']

  edge [minlen = 2]
    A->B
    A->C
    C->D
    C->E
    F->D [dir=back]
    E->G
    D->H
    E->I
  { rank = same; A; B }
  { rank = same; C }
  { rank = same; F; D; E; G}
  { rank = same; H; I }
}
")

strobe


strobe %>% 
  export_svg() %>% 
  charToRaw() %>% 
  rsvg_png("./output/prec_emp_strobe_flow.png")
