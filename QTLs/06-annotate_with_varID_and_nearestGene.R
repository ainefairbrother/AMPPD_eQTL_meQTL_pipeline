#!/usr/bin/env Rscript

# Author: Aine Fairbrother-Browne
# Date: /22
# Description: 
# 
# 
# 

#### -------------- prepare env -------------- ####

# conda activate r4-base
# screen -S 
# conda activate r4-base

#### -------------- load libs -------------- ####
require(dplyr)
require(tidyr)
require(tibble)
require(vroom)

#### -------------- import functions -------------- ####

source("~/projects/amppd_analysis/pipeline/QTLs/04-annotate_with_variantID_from_bim.R")
source("~/projects/amppd_analysis/pipeline/QTLs/05-query-api-genetics-opentargets-with-helper-funs.R")

#### -------------- implement functions -------------- ####

vroom::vroom("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/aggregated_tables/matrixeqtl_res_aggregated_p1e-03_filter.csv") %>% 
  annotate.CHR.BP.table.with.variantid(target.table=., ref.alt=F, diff.cohort="PP_PD") %>%
  apply.opentarget.query.to.variantid.column.DF(target.table=., mc.cores=20) %>% 
  vroom::vroom_write(x=., file="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/aggregated_tables/matrixeqtl_res_aggregated_p1e-03_filter_annotated.csv")