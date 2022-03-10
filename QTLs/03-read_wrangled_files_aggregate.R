#!/usr/bin/env Rscript

# Author: Aine Fairbrother-Browne
# Date: 18/22
# Description: 
# This script is run after the wrangling script 01-read_matrixeqtl_pheno_output_wrangle_write.R
# The files are further wrangled, BF is removed, snp.chr col is reinstated (to replace NA with the correct character chromosome i.e. X/Y)
# This script runs a post-wrangling for the individual cohort analysis and the mega (PP+PD) analysis
# 

#### -------------- prepare env -------------- ####

# screen -S 
# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/eqtls_mqtls
# Rscript 02-read_wrangled_files_aggregate.R

#### -------------- load libs -------------- ####
require(dplyr)
require(tidyr)
require(tibble)
require(vroom)
require(readr)
require(purrr)

#### -------------- define functions -------------- ####

read.wrangled.files.aggregate = function(file.dir, p.filter=0.05){
  
  vroom::vroom(file.dir, show_col_types = FALSE) %>% 
    # re-extract the chr, as it sometimes messes up and X/Y become NA
    tidyr::extract(col=snp.chr.pos, into="snp.chr", regex="chr(.+):", remove=FALSE) %>% 
    # if the col exists, remove it
    dplyr::select(-any_of(c("nearestGene.symbol", "nearestGene.id", "nearestGene.bioType"))) %>% 
    dplyr::filter(p.value<p.filter) %>%
    dplyr::filter(cohort!="BF") %>% 
    return(.)
}

read.wrangled.files.aggregate.meta = function(file.dir){
  
  # read in individual pheno meta files, and aggregate into single table 
  meta.aggregated = list.files(file.dir, pattern=".meta$", full.names=T) %>% 
    
    # apply to files in file.dir 
    lapply(X=., FUN=function(X){
      
      # extract relevant info from file name
      extracted.groups = stringi::stri_match_all(X, regex="PP_PD_meta_([A-Za-z]+)_pheno=(.+)_maf=0\\.05_phenotype=([A-Za-z]+)_model=([A-Za-z]+)\\.meta")[[1]]
      
      # read in file as table and add relevant cols parsed from the file name
      read.table(X, sep="", header=T) %>% 
        dplyr::mutate(
          cohort="PP_PD",
          diagnosis=extracted.groups[,2],
          phenotype=extracted.groups[,3],
          phenotype.category=extracted.groups[,4],
          matrixeqtl.model.run=extracted.groups[,5]
        ) %>% 
        return(.)
      
    }) %>% 
    dplyr::bind_rows(.) %>% 
    readr::write_csv(x=., file="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/aggregated_tables/matrixeqtl_meta_res_aggregated_no_p_filter.csv")
}

lapply.helper = function(fn, path, pattern, p.filter=0.05){
  
  # Run a function, fn, across files matching pattern, pattern, in directory, path.
  # This will produce whatever files fn is designed to output
  file.list = list.files(path=path, pattern=pattern, full.names=T)
  lapply(X=file.list, FUN=fn, p.filter=p.filter)
  
}

#### -------------- implement functions -------------- ####

# # 1. for single cohort analysis
# print("1. aggregating wrangled files for cohort analysis")
# lapply.helper(fn=read.wrangled.files.aggregate, path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/", pattern="_wrangled.csv") %>%
#   .[lapply(., length) > 0] %>%
#   .[purrr::map(., ~dim(.)[1]) > 0] %>%
#   dplyr::bind_rows() %>%
#   readr::write_csv(x=., file="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/aggregated_tables/matrixeqtl_res_aggregated_no_p_filter.csv")

# 2. for mega analysis
# print("2. aggregating wrangled files for mega analysis")
# lapply.helper(fn=read.wrangled.files.aggregate, path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/", pattern="_wrangled.csv") %>%
#   .[lapply(., length) > 0] %>%
#   .[purrr::map(., ~dim(.)[1]) > 0] %>%
#   dplyr::bind_rows() %>%
#   readr::write_csv(x=., file="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/aggregated_tables/matrixeqtl_mega_res_aggregated_no_p_filter.csv")

print("2. aggregating wrangled files for mega analysis")
lapply.helper(fn=read.wrangled.files.aggregate, path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/", pattern="_wrangled.csv") %>%
  .[lapply(., length) > 0] %>%
  .[purrr::map(., ~dim(.)[1]) > 0] %>%
  dplyr::bind_rows() %>%
  readr::write_csv(x=., file="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/aggregated_tables/matrixeqtl_mega_res_aggregated_p1e-03_filter.csv")

# # 3. for meta analysis
# print("3. aggregating wrangled files for meta analysis")
# read.wrangled.files.aggregate.meta(file.dir="/home/abrowne/projects/amppd_analysis/data/PLINK_output/")

# # 4. for timepoint analysis
# print("1. aggregating wrangled files for timepoint analysis")
# lapply.helper(fn=read.wrangled.files.aggregate, path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/", pattern="_wrangled.csv") %>%
#   .[lapply(., length) > 0] %>%
#   .[purrr::map(., ~dim(.)[1]) > 0] %>%
#   dplyr::bind_rows() %>%
#   readr::write_csv(x=., file="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/aggregated_tables/matrixeqtl_res_aggregated_no_p_filter.csv")




