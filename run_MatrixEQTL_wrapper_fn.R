#!/usr/bin/env Rscript

# use env with: conda activate r4-base

########################### Import function ###########################

source("/home/abrowne/projects/amppd_analysis/pipeline/MatrixEQTL_wrapper_fn.R")
setwd("/home/abrowne/projects/amppd_analysis/data/")

##################### 1. Run for expression data ######################

# ## run modelLINEAR ##
# 
# #for(cohort in c("PD", "PP", "BF")){
# for(cohort in c("PP")){
# 
#   for(diag in c("Case", "Control")){
# 
#     print(paste("Running modelLINEAR on expression data:", cohort, diag))
# 
#     # define lower case diag
#     if(diag=="Case"){diag.lower="case"}
#     if(diag=="Control"){diag.lower="control"}
# 
#     run.matrixeqtl(
#       SNP_file_name = paste0("./MatrixEQTL_input/genomics/celltype_PP_samples_filtered/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
#       # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
#       pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/celltype_PP_samples_filtered/",cohort,"_",diag,"_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
#       covariates_file_name = paste0("./MatrixEQTL_input/covariates/celltype_plus_standard_correction_no_PEER/",cohort,"_",diag,"_maf0.05_cov_table_for_MatrixEQTL.csv"),
#       # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
#       output_file_name = paste0("./MatrixEQTL_output/celltype_plus_standard_correction_no_PEER/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEAR_MatrixEQTL.csv"),
#       which_model="modelLINEAR",
#       which_pheno="expression"
#     )
#   }
# }

## run modelLINEAR_CROSS ##

# testing for the significance of the interaction between the genotype and the last covariate

#for(cohort in c("PD", "PP", "BF")){

# # std + celltype - peer
# for(cohort in c("PP")){
# 
#   print(paste("Running modelLINEAR_CROSS on expression data:", cohort))
# 
#   run.matrixeqtl(
#     SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
#     # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
#     pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
#     covariates_file_name = paste0("./MatrixEQTL_input/covariates/celltype_plus_standard_correction_no_PEER//",cohort,"_cohort_maf0.05_cov_table_for_MatrixEQTL.csv"),
#     # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
#     output_file_name = paste0("./MatrixEQTL_output/celltype_plus_standard_correction_no_PEER/",cohort,"_cohort_pheno=_maf=0.05_phenotype=expression_model=LINEARCROSS_MatrixEQTL.csv"),
#     which_model="modelLINEAR_CROSS",
#     which_pheno="expression"
#   )
# }

# std + celltype + peer
for(cohort in c("PP")){
  
  print(paste("Running modelLINEAR_CROSS on expression data:", cohort))
  
  run.matrixeqtl(
    SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("./MatrixEQTL_input/transcriptomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
    covariates_file_name = paste0("./MatrixEQTL_input/covariates/celltype_plus_standard_correction_with_PEER/",cohort,"_cohort_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("./MatrixEQTL_output/celltype_plus_standard_correction_with_PEER/",cohort,"_cohort_pheno=_maf=0.05_phenotype=expression_model=LINEARCROSS_MatrixEQTL.csv"),
    which_model="modelLINEAR_CROSS",
    which_pheno="expression"
  )
}

##################### 2. Run for methylation data ######################

# ## run modelLINEAR ##
# 
# #for(cohort in c("PD", "PP", "BF")){
# for(cohort in c("PP")){
# 
#   for(diag in c("Case", "Control")){
# 
#     print(paste("Running modelLINEAR on methylation data:", cohort, diag))
# 
#     # define lower case diag
#     if(diag=="Case"){diag.lower="case"}
#     if(diag=="Control"){diag.lower="control"}
# 
#     run.matrixeqtl(
#       SNP_file_name = paste0("./MatrixEQTL_input/genomics/celltype_PP_samples_filtered/",cohort,"_",diag.lower,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
#       # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
#       pheno_file_name = paste0("./MatrixEQTL_input/methylomics/celltype_PP_samples_filtered/",cohort,"_",diag,"_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
#       covariates_file_name = paste0("./MatrixEQTL_input/covariates/celltype_plus_standard_correction_no_PEER/",cohort,"_",diag,"_maf0.05_cov_table_for_MatrixEQTL.csv"),
#       # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
#       output_file_name = paste0("./MatrixEQTL_output/celltype_plus_standard_correction_no_PEER/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=methylation_model=LINEAR_MatrixEQTL.csv"),
#       which_model="modelLINEAR",
#       which_pheno="methylation"
#     )
#   }
# }

## run modelLINEAR_CROSS ##

# testing for the significance of the interaction between the genotype and the last covariate

#for(cohort in c("PD", "PP", "BF")){ 

# std + celltype - peer
for(cohort in c("PP")){ 
  
  print(paste("Running modelLINEAR_CROSS on methylation data:", cohort))
  
  run.matrixeqtl(
    SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("./MatrixEQTL_input/methylomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
    covariates_file_name = paste0("./MatrixEQTL_input/covariates/celltype_plus_standard_correction_no_PEER/",cohort,"_cohort_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("./MatrixEQTL_output/celltype_plus_standard_correction_no_PEER/",cohort,"_cohort_pheno=_maf=0.05_phenotype=methylation_model=LINEARCROSS_MatrixEQTL.csv"),
    which_model="modelLINEAR_CROSS",
    which_pheno="methylation"
  )
}

# std + celltype + peer
for(cohort in c("PP")){ 
  
  print(paste("Running modelLINEAR_CROSS on methylation data:", cohort))
  
  run.matrixeqtl(
    SNP_file_name = paste0("./MatrixEQTL_input/genomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv"),
    # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    pheno_file_name = paste0("./MatrixEQTL_input/methylomics/PP_filtered_for_samples_with_celltype_data/",cohort,"_cohort_pheno=_methylation_matrix_MatrixEQTL_input.csv"),
    covariates_file_name = paste0("./MatrixEQTL_input/covariates/celltype_plus_standard_correction_with_PEER/",cohort,"_cohort_maf0.05_cov_table_for_MatrixEQTL.csv"),
    # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
    output_file_name = paste0("./MatrixEQTL_output/celltype_plus_standard_correction_with_PEER/",cohort,"_cohort_pheno=_maf=0.05_phenotype=methylation_model=LINEARCROSS_MatrixEQTL.csv"),
    which_model="modelLINEAR_CROSS",
    which_pheno="methylation"
  )
}
