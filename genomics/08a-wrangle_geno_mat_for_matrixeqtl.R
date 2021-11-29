#!/usr/bin/env Rscript

for(cohort in c("BF", "PD", "PP")){

# import libs
require(data.table)
require(tidyverse)

# set file suffixes
mat.suffix="_all_chrs-maf0.05_geno_matrix.012"
pos.suffix="_all_chrs-maf0.05_geno_matrix.012.pos"
indv.suffix="_all_chrs-maf0.05_geno_matrix.012.indv"
out.suffix="_all_chrs-maf0.05_geno_matrix_012_wrangled.csv"

# define directories
base.dir="/scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files/"
sample.id.dir="/scratch/groups/hodgkinsonlab/aine/data/AMPPD_sample_lists/"

# import mat
print("Importing geno_matrix.012")
mat = data.table::fread(paste0(base.dir,cohort,mat.suffix))

# remove V1 as it's an index col
mat = mat[, -"V1"]

# check dim 
mat %>% dim() %>% print()

# import the indv (samples) and pos (genomic pos) files
print("Importing 012.indv and 012.pos")
indv=data.table::fread(paste0(base.dir,cohort,indv.suffix), header=F) %>% dplyr::pull(V1)

# import pos and make character col label vector by pasting chr and numeric position together
pos=data.table::fread(paste0(base.dir,cohort,pos.suffix), header=F) %>% dplyr::mutate(V3 = paste0(V1, "-", V2)) %>% dplyr::pull(V3)

# check lengths of indv and pos
print("Check lengths of indv and pos")
indv %>% length() %>% print()
pos %>% length() %>% print()

# transpose mat so that cols=samples rows=pos
mat = mat %>% t()

# assign indv as colnames to the geno matrix
colnames(mat) = indv

# assign pos as rownames to the geno matrix
rownames(mat) = pos

# check mat 
print("Checking rownames, colnames, dim and head")
rownames(mat)[1:5] %>% print()
colnames(mat)[1:5] %>% print()
dim(mat) %>% print()
mat[1:10,1:10] %>% print()

# import 1015 sample lists
print("Importing sample lists")

samples.to.use.1015 = readr::read_csv(paste0(sample.id.dir,cohort,"-visit_month=0-case_control_other_latest!=Other_unique_1015.csv"), show_col_types = FALSE) %>% 
  dplyr::pull(x) %>% 
  unique() %>%
  gsub("-BLM0T1", "", .) %>% 
  gsub("-SVM0_5T1", "", .)

print(paste("length of samples.to.use.1015=", length(samples.to.use.1015)))

samples.to.use.1015.case = readr::read_csv(paste0(sample.id.dir,cohort,"-visit_month=0-case_control_other_latest==Case_unique_1015.csv"), show_col_types = FALSE) %>% 
  dplyr::pull(x) %>% 
  unique() %>%
  gsub("-BLM0T1", "", .) %>% 
  gsub("-SVM0_5T1", "", .)

print(paste("length of samples.to.use.1015.case=", length(samples.to.use.1015.case)))

samples.to.use.1015.control = readr::read_csv(paste0(sample.id.dir,cohort,"-visit_month=0-case_control_other_latest==Control_unique_1015.csv"), show_col_types = FALSE) %>% 
  dplyr::pull(x) %>% 
  unique() %>%
  gsub("-BLM0T1", "", .) %>% 
  gsub("-SVM0_5T1", "", .)

print(paste("length of samples.to.use.1015.control=", length(samples.to.use.1015.control)))

# filter and write out
print("Writing data...")

mat %>% .[, colnames(.) %in% samples.to.use.1015] %>% as.data.frame() %>% tibble::rownames_to_column("gene_id") %>% head() %>% print()

fwrite(x=mat %>% .[, colnames(.) %in% samples.to.use.1015] %>% as.data.frame() %>% tibble::rownames_to_column("gene_id"), file=paste0(base.dir,cohort,"_cohort_",out.suffix))

fwrite(x=mat %>% .[, colnames(.) %in% samples.to.use.1015.case] %>% as.data.frame() %>% tibble::rownames_to_column("gene_id"), file=paste0(base.dir,cohort,"_Case_",out.suffix))

fwrite(x=mat %>% .[, colnames(.) %in% samples.to.use.1015.control] %>%  as.data.frame()  %>% tibble::rownames_to_column("gene_id"), file=paste0(base.dir,cohort,"_Control_",out.suffix))

rm(list = ls())

}


