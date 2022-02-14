#!/usr/bin/env Rscript

# Author: Aine Fairbrother-Browne
# Date: 02/22
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

#### -------------- define functions -------------- ####

# Read in and filter for samples in AMP-PD release `2019_v1release_1015`  




pp.samples.to.use = read_csv("./data/AMPPD_sample_lists/PP-visit_month=0-case_control_other_latest!=Other_unique_1015.csv", show_col_types = FALSE) %>% 
  dplyr::pull(x)

pp = read.csv("./data/mt_nuc_corr_pipeline_output/PP_log10_mediannorm_TPM.csv", header=F) 
colnames(pp) = c("gene_id", pp[1,] %>% as.character()) %>% .[-length(.)]

pp = pp %>% 
  tidyr::drop_na() %>% 
  tibble::column_to_rownames("gene_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate_if(is.character, as.numeric) %>% 
  `rownames<-`(gsub("\\.", "-", rownames(.))) %>% 
  .[rownames(.) %in% pp.samples.to.use, ]

# Add gene symbols, calculate quartiles  

pp.long = pp %>% 
  tibble::rownames_to_column("sample_id") %>% 
  tidyr::pivot_longer(2:ncol(.)) %>% 
  dplyr::mutate(GENEID=gsub("\\..*", "",name)) 

pp.long.quartiles = pp.long %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::mutate(q1=quantile(value, 0.25, na.rm=T),
                q3=quantile(value, 0.75, na.rm=T)) %>% 
  dplyr::mutate(iqr=q3-q1) %>% 
  dplyr::mutate(outlier.low=q1-(iqr*3), 
                outlier.high=q3+(iqr*3))

# Mask outliers and write out 

pp.long.quartiles %>% 
  dplyr::mutate(
    is.outlier = case_when(
      value > outlier.high ~ TRUE,
      value < outlier.low ~ TRUE,
      (value <= outlier.high) & (value >= outlier.low) ~ FALSE)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(value = ifelse(is.outlier==TRUE, NA, value)) %>% 
  dplyr::select(sample_id, name, value) %>%
  tidyr::pivot_wider(id_cols=sample_id, names_from=name, values_from=value) %>%
  write_csv(x=., file="./data/mt_nuc_corr_pipeline_output/PP_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr.csv")








#### -------------- implement functions -------------- ####
