#!/usr/bin/env Rscript

########################### Purpose #######################################

# This wrapper script runs MatrixEQTL for the individual cohort analysis

########################### Set up environment ############################

# screen -S run-tp-avg
# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/
# Rscript run_MatrixEQTL_wrapper_timepoint_avg.R >& MatrixEQTL_log.txt
 
########################### Import functions ###########################

source("/home/abrowne/projects/amppd_analysis/pipeline/MatrixEQTL_wrapper_fn.R")
setwd("/home/abrowne/projects/amppd_analysis/data/")

########################### set up script timer ###########################

start_time = Sys.time()

##################### 1. Run for expression data ######################

## run modelLINEAR ##
for(cohort in c("PP", "PD")){
  for(diag in c("Case", "Control")){

    print(paste("Running modelLINEAR on expression data:", cohort, diag))

    # define lower case diag
    if(diag=="Case"){diag.lower="case"}
    if(diag=="Control"){diag.lower="control"}

    run.matrixeqtl(
      SNP_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/all_timepoints/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
      # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
      pheno_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/transcriptomics/all_timepoints/",cohort,"_",diag,"_timepoint_averaged_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
      covariates_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/all_timepoints/",cohort,"_",diag,"_averaged_timepoint_maf0.05_cov_table_for_MatrixEQTL.csv"),
      # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
      output_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEAR_MatrixEQTL.csv"),
      which_model="modelLINEAR",
      which_pheno="expression"
    )
  }
}

## run modelLINEAR_CROSS ##
# testing for the significance of the interaction between the genotype and the last covariate
for(cohort in c("PD", "PP")){
  
  print(paste("Running modelLINEAR_CROSS on expression data:", cohort))
  diag="cohort"

  run.matrixeqtl(
    SNP_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/all_timepoints/",cohort,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/transcriptomics/all_timepoints/",cohort,"_",diag,"_timepoint_averaged_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
    covariates_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/all_timepoints/",cohort,"_",diag,"_averaged_timepoint_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEARCROSS_MatrixEQTL.csv"),
    which_model="modelLINEAR_CROSS",
    which_pheno="expression"
  )
}

##################### 2. Run for methylation data ######################

## run modelLINEAR ##
for(cohort in c("PP", "PD")){
  for(diag in c("Control", "Case")){
    
    print(paste("Running modelLINEAR on methylation data:", cohort, diag))
    
    # define lower case diag
    if(diag=="Case"){diag.lower="case"}
    if(diag=="Control"){diag.lower="control"}
    
    run.matrixeqtl(
      SNP_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/all_timepoints/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
      # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
      pheno_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/methylomics/all_timepoints/",cohort,"_",diag,"_all_timepoints_avg_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
      covariates_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/all_timepoints/",cohort,"_",diag,"_averaged_timepoint_maf0.05_cov_table_for_MatrixEQTL.csv"),
      # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
      output_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=methylation_model=LINEAR_MatrixEQTL.csv"),
      which_model="modelLINEAR",
      which_pheno="methylation"
    )
  }
}

## run modelLINEAR_CROSS ##
# testing for the significance of the interaction between the genotype and the last covariate
for(cohort in c("PP", "PD")){
  
  print(paste("Running modelLINEAR_CROSS on methylation data:", cohort))
  diag="cohort"
  
  run.matrixeqtl(
    SNP_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/all_timepoints/",cohort,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/methylomics/all_timepoints/",cohort,"_cohort_all_timepoints_avg_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
    covariates_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/all_timepoints/",cohort,"_cohort_averaged_timepoint_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=methylation_model=LINEARCROSS_MatrixEQTL.csv"),
    which_model="modelLINEAR_CROSS",
    which_pheno="methylation"
  )
}

########################### get time taken for script to run ###########################

end_time = Sys.time()
time_taken = end_time - start_time
print(time_taken)

########################### end of script ###########################
