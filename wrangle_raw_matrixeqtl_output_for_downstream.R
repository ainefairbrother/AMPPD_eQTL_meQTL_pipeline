#!/usr/bin/env Rscript

# wrangling raw matrixeqtl output, preparing it for downstream 

#### -------------- prepare env -------------- ####

# conda activate r4-base

# conda install -c conda-forge r-stringi
# conda install -c conda-forge r-tzdb
# conda install -c r r-tidyverse

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
    extracted.groups = stringi::stri_match_all(X, regex="^([A-Z]{2})_([A-Za-z]+)_pheno=(.+)_maf=(.+)_phenotype=([A-Za-z]+)_model=([A-Za-z]+)_MatrixEQTL.csv")[[1]]
    
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

read.annot.file.apply.filter = function(file.path, p.cutoff=1e-03){
  
  d = vroom::vroom(file.path) 
  
  if( !("nearestGene.id" %in% colnames(d))) {
    d = d %>% 
      dplyr::mutate(nearestGene.symbol=NA, 
                    nearestGene.id=NA,
                    nearestGene.bioType=NA)
  }
  d %>%
    dplyr::filter(p.value<p.cutoff) %>%
    # make snp col a number
    dplyr::mutate(snp.chr=as.numeric(snp.chr)) %>% 
    return(.)
}
lapply.helper = function(fn, path, pattern){
  
  # Run a function, fn, across files matching pattern, pattern, in directory, path. 
  # This will produce whatever files fn is designed to output
  file.list = list.files(path=path, pattern=pattern, full.names=T)
  lapply(X=file.list, FUN=fn)
  
}

#### -------------- implement functions -------------- ####

## implement read.matrixeqtl.pheno.output.wrangle.write
# read.matrixeqtl.pheno.output.wrangle.write(path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/celltype_plus_standard_correction_no_PEER/", 
#                                            pattern="MatrixEQTL.csv")

# read.matrixeqtl.pheno.output.wrangle.write(path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/celltype_plus_standard_correction_with_PEER/", 
#                                            pattern="MatrixEQTL.csv")

## implement read.annot.file.apply.filter on wrangled raw files, collate into aggregated table and write out
lapply.helper(fn=read.annot.file.apply.filter, path="./data/MatrixEQTL_output", pattern="_wrangled.csv", p.cutoff=1e-03) %>% 
  # then filter list of files and bind them together
  .[lapply(., length) > 0] %>% 
  .[map(., ~dim(.)[1]) > 0] %>% 
  dplyr::bind_rows() %>% 
  # remove extra cols 
  dplyr::select(-(ends_with(".x"))) %>% 
  dplyr::select(-(ends_with(".y"))) %>% 
  dplyr::select(-contains("nearestGene")) %>% 
  readr::write_csv(x=., file="./data/MatrixEQTL_output/aggregated_tables/matrixeqtl_res_aggregated_p1e-03_filter.csv")









