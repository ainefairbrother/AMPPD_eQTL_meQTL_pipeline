#!/usr/bin/env Rscript

# Author: Aine Fairbrother-Browne
# Date: 02/22
# Description: 
# 
# 
# 

#### -------------- prepare env -------------- ####

# screen -S calc-outliers
# conda activate r4-base
# /home/abrowne/projects/amppd_analysis/pipeline/transcriptomics
# 08-mask_outliers_q1-3iqr_q3+3iqr_for_log10mednorm_files.R

#### -------------- load libs -------------- ####
require(dplyr)
require(tidyr)
require(tibble)
require(readr)
require(vroom)
require(tibble)
require(IRanges)

setwd("/home/abrowne/projects/amppd_analysis/")

#### -------------- define functions -------------- ####

# Read in and filter for samples in AMP-PD release `2019_v1release_1015`  

run.mask.outliers = function(path.to.sample.list, path.to.log10mednorm, path.out){
  
  samples.to.use = readr::read_csv(path.to.sample.list, show_col_types = FALSE) %>% 
    tidyr::separate(col=x, into=c("cohort", "participant", "visit"), remove=F) %>% 
    dplyr::mutate(participant_id=paste0(cohort,"-",participant)) %>% 
    dplyr::pull(participant_id)
  
  # test
  #vroom::vroom(path.to.log10mednorm) %>% dplyr::select(-gene_id) %>% anyNA() %>% print()
  vroom::vroom(path.to.log10mednorm) %>% 
    tidyr::drop_na() %>% 
    tibble::column_to_rownames("gene_id") %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::mutate_if(is.character, as.numeric) %>% 
    .[rownames(.) %in% samples.to.use, ] %>% 
    
    # calculate quartiles  
    tibble::rownames_to_column("sample_id") %>% 
    tidyr::pivot_longer(2:ncol(.)) %>% 
    dplyr::group_by(sample_id) %>% 
    dplyr::mutate(q1=quantile(value, 0.25, na.rm=T),
                  q3=quantile(value, 0.75, na.rm=T)) %>%
    dplyr::mutate(iqr=q3-q1) %>% 
    dplyr::mutate(outlier.low=q1-(iqr*3), 
                  outlier.high=q3+(iqr*3)) %>% 
    
    # Mask outliers and write out 
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
    # test 
    # dplyr::select(-sample_id) %>% anyNA() %>% print()
    readr::write_csv(x=., file=path.out)
}

#### -------------- implement functions -------------- ####

for(cohort in c("PP", "PD")){
  for(diag in c("case", "ctrl")){
    
    if(diag=="case"){diag.upper="Case"}
    if(diag=="ctrl"){diag.upper="Control"}
    
    run.mask.outliers(
      path.to.sample.list=paste0("./data/AMPPD_sample_lists/",cohort,"-visit_month=ALL-case_control_other_latest==",diag.upper,"_unique_1015.csv"),
      path.to.log10mednorm=paste0("./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/",cohort,"_",diag,"_timepoint_averaged_log10_mediannorm_TPM.csv"),
      path.out=paste0("./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/",cohort,"_",diag.upper,"_timepoint_averaged_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr.csv")
    )
  }
  for(diag in c("cohort")){
    run.mask.outliers(
      path.to.sample.list=paste0("./data/AMPPD_sample_lists/",cohort,"-visit_month=ALL-case_control_other_latest!=Other_unique_1015.csv"),
      path.to.log10mednorm=paste0("./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/",cohort,"_",diag,"_timepoint_averaged_log10_mediannorm_TPM.csv"),
      path.out=paste0("./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/",cohort,"_",diag,"_timepoint_averaged_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr.csv")
    )
  }
}








