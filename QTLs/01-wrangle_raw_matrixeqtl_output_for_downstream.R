#!/usr/bin/env Rscript

# wrangling raw matrixeqtl output, preparing it for downstream 

#### -------------- prepare env -------------- ####

# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/

# the following had to be installed:
# conda install -c conda-forge r-stringi
# conda install -c conda-forge r-tzdb
# conda install -c r r-tidyversecd

#### -------------- load libs -------------- ####
library(dplyr)
library(tidyr)
library(tibble)
library(stringi)
library(vroom)
library(parallel)
library(grDevices)
library(stats)

#### -------------- define functions -------------- ####

# 1. Basic wrangling. Raw -> wrangled.
# MatrixEQTL.csv -> _wrangled.csv
# extract file name information into variables, 
# generate relevant variables like chr.pos and snp.chr.pos
# assign datatypes
# calculate beta_se
read.matrixeqtl.pheno.output.wrangle.write = function(path, pattern){
  
  setwd(path)
  
  # get all the files to wrangle
  file.list = list.files(path=path, pattern=pattern, full.names=F)
  
  # define function to apply to each matrixeqtl output file in file.list
  wrangle.write.file=function(X){
    
    # for stop-start running - don't generate file if it already exists
    #if(file.exists(gsub("MatrixEQTL", "wrangled", X))==FALSE){
    
    # extract relevant info from file name
    extracted.groups = stringi::stri_match_all(X, regex="(.+)_([A-Za-z]+)_pheno=(.+)_maf=(.+)_phenotype=([A-Za-z]+)_model=([A-Za-z]+)_MatrixEQTL.csv")[[1]]
    
    # read in and wrangle - add file name info as columns 
    vroom::vroom(file=X, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
      dplyr::mutate(
        beta=as.numeric(beta),
        `t-stat`=as.numeric(`t-stat`),
        `p-value`=as.numeric(`p-value`),
        `FDR`=as.numeric(`FDR`),
        cohort=rep(extracted.groups[,2], nrow(.)),
        diagnosis=extracted.groups[,3],
        phenotype=extracted.groups[,4],
        nuc.snp.maf.filter=extracted.groups[,5],
        phenotype.category=extracted.groups[,6],
        matrixeqtl.model.run=extracted.groups[,7]) %>% 
      
      # do some additional tidying of cols
      tidyr::separate(col=SNP, into=c("chr", "start"), sep="-") %>% 
      dplyr::select(-gene) %>% 
      dplyr::mutate(chr.pos = paste0(chr, ":", start)) %>% 
      dplyr::relocate(chr.pos, chr, start, phenotype) %>% 
      dplyr::rename(snp.chr.pos=chr.pos, snp.chr=chr, snp.start=start) %>% 
      dplyr::rename(p.value = `p-value`, t.stat = `t-stat`) %>%  
      dplyr::mutate(snp.chr=gsub("chr", "", snp.chr)) %>% 
      dplyr::mutate(snp.start=as.numeric(snp.start)) %>% 
      dplyr::arrange(p.value) %>% 
      # calculate the standard error of the beta - useful for downstream meta-analyses
      dplyr::mutate(beta_se = beta/t.stat) %>% 
      # write out
      vroom::vroom_write(x=., file=gsub("MatrixEQTL", "wrangled", X), delim=",")
  }
  parallel::mclapply(file.list, wrangle.write.file, mc.cores=4)
}

#### -------------- implement functions -------------- ####

## implement read.matrixeqtl.pheno.output.wrangle.write
# read.matrixeqtl.pheno.output.wrangle.write(path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/celltype_plus_standard_correction_no_PEER/", 
#                                            pattern="MatrixEQTL.csv")

# read.matrixeqtl.pheno.output.wrangle.write(path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/celltype_plus_standard_correction_with_PEER/", 
#                                            pattern="MatrixEQTL.csv")

# # run for mega
read.matrixeqtl.pheno.output.wrangle.write(path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/",
                                           pattern="MatrixEQTL.csv")

# # run for timepoint
# read.matrixeqtl.pheno.output.wrangle.write(path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/",
#                                            pattern="MatrixEQTL.csv")



