#!/usr/bin/env bash

# Author: Aine Fairbrother-Browne
# Date: 02/22
# Description: 
# this script is to wrangle the matrixeqtl input data to run a mega-analysis
# PD and PP cohorts will be combined to prepare data for a combined-cohort mega-analysis
# so this script will combine the transcriptomic, methylomic, genomic and covariate data of the two cohorts 

#### -------------- prepare env -------------- ####

# screen -S wrangle-mega-timepoint
# conda activate r4-base
# cd ~/projects/amppd_analysis/pipeline
# Rscript wrangle wrangle_MatrixEQTL_input_data_for_PP_PD_timepoint_megaanalysis.R

#### -------------- load libs -------------- ####
require(dplyr)
require(tidyr)
require(tibble)
require(vroom)

#### -------------- define functions -------------- ####

# these funtions take in wrangled matrixeqtl input files, merge the PP and PD cohorts
# PP_PD files are output into /in/path/PP_PD_mega_analysis/

wrangle.covariate.data.for.megaanalysis = function(infile.path){
  
  outfile.path=paste0(gsub("all_timepoints/", "", infile.path),"PP_PD_timepoint_mega_analysis/")
  suffix.list = list.files(path=infile.path, pattern="_averaged_timepoint_maf0.05_cov_table_for_MatrixEQTL.csv") %>% 
    gsub(pattern="^\\w{2}_", "", .) %>% 
    unique()
  
  print(suffix.list)
  
  for(suffix in suffix.list){
    dplyr::bind_cols(
      vroom::vroom(paste0(infile.path, "PP_", suffix), show_col_types = FALSE),
      vroom::vroom(paste0(infile.path, "PD_", suffix), show_col_types = FALSE) %>%
        dplyr::select(-variable)
    ) %>% 
      vroom::vroom_write(., paste0(outfile.path,"PP_PD_",suffix), delim = ",")
  }
}
wrangle.genomic.data.for.megaanalysis = function(infile.path){
  
  outfile.path=paste0(gsub("all_timepoints/", "", infile.path),"PP_PD_timepoint_mega_analysis/")
  suffix.list = list.files(path=infile.path, pattern="_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv") %>% 
    gsub(pattern="^\\w{2}_", "", .) %>% 
    unique()
  
  print(suffix.list)
  
  for(suffix in suffix.list){
    
    pos.in.both=intersect(vroom::vroom(paste0(infile.path,"PP_",suffix), show_col_types = FALSE, col_select=c(pos))$pos, 
                          vroom::vroom(paste0(infile.path,"PD_",suffix), show_col_types = FALSE, col_select=c(pos))$pos
    )
    
    dplyr::bind_cols(
      
      vroom::vroom(paste0(infile.path,"PP_",suffix), show_col_types = FALSE) %>% 
        dplyr::filter(pos %in% pos.in.both) %>% 
        dplyr::arrange(match(pos, pos.in.both)),
      
      vroom::vroom(paste0(infile.path,"PD_",suffix), show_col_types = FALSE) %>% 
        dplyr::filter(pos %in% pos.in.both) %>% 
        dplyr::arrange(match(pos, pos.in.both)) %>% 
        dplyr::select(-pos)
      
    ) %>% 
      vroom::vroom_write(., paste0(outfile.path,"PP_PD_",suffix), delim = ",")
  }
}
wrangle.transcriptomic.data.for.megaanalysis = function(infile.path){
  
  outfile.path=paste0(gsub("all_timepoints/", "", infile.path),"PP_PD_timepoint_mega_analysis/")
  suffix.list = list.files(path=infile.path, pattern="P*_pheno=ENSG\\d*_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3\\+3iqr_MT_ONLY.csv") %>% 
                                                  # PP_Control_timepoint_averaged_pheno=ENSG00000212907_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv
    gsub(pattern="^\\w{2}_", "", .) %>% 
    unique()
  
  suffix.list = suffix.list[!grepl("ENSG00000228253", suffix.list)]
  print(suffix.list)
  
  for(suffix in suffix.list){
    
    dplyr::bind_cols(
      vroom::vroom(paste0(infile.path, "PP_", suffix), show_col_types = FALSE),
      vroom::vroom(paste0(infile.path, "PD_", suffix), show_col_types = FALSE) %>%
        dplyr::select(-sample_id)
    ) %>% 
      vroom::vroom_write(., paste0(outfile.path,"PP_PD_",suffix), delim = ",")
  }
}
wrangle.methylomic.data.for.megaanalysis = function(infile.path){
  
  outfile.path=paste0(gsub("all_timepoints/", "", infile.path),"PP_PD_timepoint_mega_analysis/")
  suffix.list = list.files(path=infile.path, pattern="P*_pheno=\\d*_full_methylation_matrix_MatrixEQTL_input.csv") %>% 
    gsub(pattern="^\\w{2}_", "", .) %>% 
    unique()
  
  print(suffix.list)
  
  for(suffix in suffix.list){
    
    dplyr::bind_cols(
      vroom::vroom(paste0(infile.path, "PP_", suffix), show_col_types = FALSE),
      vroom::vroom(paste0(infile.path, "PD_", suffix), show_col_types = FALSE) %>%
        dplyr::select(-sample_id)
    ) %>% 
      vroom::vroom_write(., paste0(outfile.path,"PP_PD_",suffix), delim = ",")
  }
}

#### -------------- implement functions -------------- ####

print("wrangling covariate data")
wrangle.covariate.data.for.megaanalysis(infile.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/all_timepoints/")
print("wrangling transcriptomic data")
wrangle.transcriptomic.data.for.megaanalysis(infile.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/transcriptomics/all_timepoints/")
print("wrangling methylomic data")
wrangle.methylomic.data.for.megaanalysis(infile.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/methylomics/all_timepoints/")
print("wrangling genomic data")
wrangle.genomic.data.for.megaanalysis(infile.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/all_timepoints/")

