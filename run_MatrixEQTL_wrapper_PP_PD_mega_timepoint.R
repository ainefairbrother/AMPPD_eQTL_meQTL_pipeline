#!/usr/bin/env Rscript

########################### Purpose #######################################

# This wrapper script runs MatrixEQTL for the mega timepoint analysis

########################### Set up environment ############################

# screen -S run-mega-tp
# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/
# Rscript run_MatrixEQTL_wrapper_PP_PD_mega_timepoint.R

########################### Import functions ###########################

source("/home/abrowne/projects/amppd_analysis/pipeline/MatrixEQTL_wrapper_fn.R")
setwd("/home/abrowne/projects/amppd_analysis/data/")

# define cohort label 
cohort="PP_PD"

########################### set up script timer ###########################

start_time = Sys.time()

##################### 1. Run for expression data ######################

## run modelLINEAR ##
for(diag in c("Case", "Control")){

  print(paste("Running modelLINEAR on expression data:", cohort, diag))

  # define lower case diag
  if(diag=="Case"){diag.lower="case"}
  if(diag=="Control"){diag.lower="control"}

  run.matrixeqtl(
    SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_PD_mega_analysis/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/PP_PD_timepoint_mega_analysis/",cohort,"_",diag,"_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
    covariates_file_name = paste0("./MatrixEQTL_input/covariates/PP_PD_timepoint_mega_analysis/",cohort,"_",diag,"_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("./MatrixEQTL_output/PP_PD_timepoint_mega_analysis/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEAR_MatrixEQTL.csv"),
    which_model="modelLINEAR",
    which_pheno="expression"
  )
}

## run modelLINEAR_CROSS ##
# testing for the significance of the interaction between the genotype and the last covariate
print(paste("Running modelLINEAR_CROSS on expression data:", cohort))
run.matrixeqtl(
  SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_PD_mega_analysis/",cohort,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
  # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/PP_PD_timepoint_mega_analysis/",cohort,"_cohort_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
  covariates_file_name = paste0("./MatrixEQTL_input/covariates/PP_PD_timepoint_mega_analysis/",cohort,"_cohort_maf0.05_cov_table_for_MatrixEQTL.csv"),
  # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  output_file_name = paste0("./MatrixEQTL_output/PP_PD_timepoint_mega_analysis/",cohort,"_cohort_pheno=_maf=0.05_phenotype=expression_model=LINEARCROSS_MatrixEQTL.csv"),
  which_model="modelLINEAR_CROSS",
  which_pheno="expression"
)


##################### 2. Run for methylation data ######################

## run modelLINEAR ##
for(diag in c("Case", "Control")){

  print(paste("Running modelLINEAR on methylation data:", cohort, diag))

  # define lower case diag
  if(diag=="Case"){diag.lower="case"}
  if(diag=="Control"){diag.lower="control"}

  run.matrixeqtl(
    SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_PD_mega_analysis/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("./MatrixEQTL_input/methylomics/PP_PD_timepoint_mega_analysis/",cohort,"_",diag,"_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
    covariates_file_name = paste0("./MatrixEQTL_input/covariates/PP_PD_timepoint_mega_analysis/",cohort,"_",diag,"_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("./MatrixEQTL_output/PP_PD_timepoint_mega_analysis/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=methylation_model=LINEAR_MatrixEQTL.csv"),
    which_model="modelLINEAR",
    which_pheno="methylation"
  )
}

## run modelLINEAR_CROSS ##
# testing for the significance of the interaction between the genotype and the last covariate
print(paste("Running modelLINEAR_CROSS on methylation data:", cohort))
run.matrixeqtl(
  SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_PD_mega_analysis/",cohort,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"), # use same as individual cohort mega analysis 
  # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  pheno_file_name = paste0("./MatrixEQTL_input/methylomics/PP_PD_timepoint_mega_analysis/",cohort,"_cohort_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
  covariates_file_name = paste0("./MatrixEQTL_input/covariates/PP_PD_timepoint_mega_analysis/",cohort,"_cohort_maf0.05_cov_table_for_MatrixEQTL.csv"),
  # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  output_file_name = paste0("./MatrixEQTL_output/PP_PD_timepoint_mega_analysis/",cohort,"_cohort_pheno=_maf=0.05_phenotype=methylation_model=LINEARCROSS_MatrixEQTL.csv"),
  which_model="modelLINEAR_CROSS",
  which_pheno="methylation"
)

########################### get time taken for script to run ###########################

end_time = Sys.time()
time_taken = end_time - start_time
print(time_taken)

########################### end of script ###########################
