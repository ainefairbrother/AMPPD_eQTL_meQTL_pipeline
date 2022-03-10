#!/usr/bin/env Rscript

########################### Purpose #######################################

# This wrapper script runs MatrixEQTL to test the transcriptomic processing pipeline 

########################### Set up environment ############################

# screen -S 
# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/
# Rscript 

########################### Import functions ###########################

source("/home/abrowne/projects/amppd_analysis/pipeline/MatrixEQTL_wrapper_fn.R")
setwd("/home/abrowne/projects/amppd_analysis/data/")


########################### set up script timer ###########################

start_time = Sys.time()

##################### 1. Run for expression data ######################

# define set params
diag = "Control"
diag.lower = "control"
cohort = "PP"

## run for PP_ctrl_log10_norm_TPM_MTonly_masked_outliers_wrangled.csv ###

## run modelLINEAR ##
print(paste("Running modelLINEAR on expression data:", cohort, diag))

run.matrixeqtl(
  SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
  # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/PP_Control_testing_transcriptomic_processing/PP_ctrl_log10_norm_TPM_MTonly_masked_outliers_wrangled.csv",cohort,"_",diag,""),
  covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_",diag,"_maf0.05_cov_table_for_MatrixEQTL.csv"),
  # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  output_file_name = paste0("./MatrixEQTL_output/PP_Control_testing_transcriptomic_processing/log10_TPM_MT_genes_masked/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEAR_MatrixEQTL.csv"),
  which_model="modelLINEAR",
  which_pheno="expression"
)

## run for PP_ctrl_log10_norm_TPM_MTonly_masked_outliers_wrangled.csv ###

## run modelLINEAR ##
print(paste("Running modelLINEAR on expression data:", cohort, diag))

run.matrixeqtl(
  SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
  # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/PP_Control_testing_transcriptomic_processing/PP_ctrl_log10_norm_TPM_MTonly_mednorm_masked_outliers_wrangled.csv",cohort,"_",diag,""),
  covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_",diag,"_maf0.05_cov_table_for_MatrixEQTL.csv"),
  # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  output_file_name = paste0("./MatrixEQTL_output/PP_Control_testing_transcriptomic_processing/log10_TPM_MT_genes_mednorm_masked/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEAR_MatrixEQTL.csv"),
  which_model="modelLINEAR",
  which_pheno="expression"
)


########################### get time taken for script to run ###########################

end_time = Sys.time()
time_taken = end_time - start_time
print(time_taken)

########################### end of script ###########################
