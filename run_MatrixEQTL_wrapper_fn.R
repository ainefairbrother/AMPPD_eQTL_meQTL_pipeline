########################### Import function ###########################

source("/home/abrowne/projects/amppd_analysis/pipeline/MatrixEQTL_wrapper_fn.R")

##################### 1. Run for expression data ######################

expression.phenos = c("ENSG00000198888", "ENSG00000198763", "ENSG00000198840",
                      "ENSG00000212907", "ENSG00000198886", "ENSG00000198786", 
                      "ENSG00000198695","ENSG00000198712", "ENSG00000198938", 
                      "ENSG00000228253", "ENSG00000198727", "ENSG00000198804", 
                      "ENSG00000198899", "ENSG00000211459","ENSG00000210082")

## run modelLINEAR ##

for(cohort in c("PD", "PP", "BF")){
  
  for(diag in c("Case", "Control")){
    
    # define lower case diag
    if(diag=="Case"){diag.lower="case"}
    if(diag=="Control"){diag.lower="control"}
    
    for(pheno in expression.phenos){
      
      run.matrixeqtl(
        SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_",diag.lower,"_all_chrs_geno_matrix_012_wrangled-8b.csv"),
        expression_file_name = paste0("./MatrixEQTL_input/transcriptomics/",cohort,"_",diag,"_pheno=",pheno,"_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
        covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_",diag,"_cov_table_for_MatrixEQTL.csv"),
        output_file_name = paste0("./MatrixEQTL_output/",cohort,"_",diag,"_pheno=",pheno,"_MatrixEQTL_pheno=expression_model=LINEAR.csv"),
        which_model="modelLINEAR"
      )
    }
  }
}

## run modelLINEAR_CROSS ##

# testing for the significance of the interaction between the genotype and the last covariate

for(cohort in c("PD", "PP", "BF")){
  
  for(pheno in expression.phenos){
    
    run.matrixeqtl(
      SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_all_chrs_geno_matrix_012_wrangled-8b.csv"),
      expression_file_name = paste0("./MatrixEQTL_input/transcriptomics/",cohort,"_cohort_pheno=",pheno"_cohort_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"),
      covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_cohort_cov_table_for_MatrixEQTL.csv"),
      output_file_name = paste0("./MatrixEQTL_output/",cohort,"_pheno=",pheno,"_cohort_MatrixEQTL_pheno=expression_model=LINEARCROSS.csv"),
      which_model="modelLINEAR_CROSS"
    )
  }
}

##################### 2. Run for methylation data ######################

# define all phenotypes in methylation data
methylation.phenos = c("3238_full", "4271_full", "4392_full", "5520_full", 
                       "5647_full", "5721_full", "5818_full", "5883_full", 
                       "7526_full", "8303_full", "9999_full", "10413_full", 
                       "12146_full", "12274_full", "14734_full", "15896_full", 
                       "15948_full", "2617_full", "13710_full")

## run modelLINEAR ##

for(cohort in c("PD", "PP", "BF")){
  
  for(diag in c("Case", "Control")){
    
    # define lower case diag
    if(diag=="Case"){diag.lower="case"}
    if(diag=="Control"){diag.lower="control"}
    
    for(pheno in methylation.phenos){
      
      run.matrixeqtl(
        SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_",diag.lower,"_all_chrs_geno_matrix_012_wrangled-8b.csv"),
        expression_file_name = paste0("./MatrixEQTL_input/methylomics/",cohort,"_",diag,"_pheno=",pheno,"_methylation_matrix_MatrixEQTL_input.csv"),
        covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_",diag,"_cov_table_for_MatrixEQTL.csv"),
        output_file_name = paste0("./MatrixEQTL_output/",cohort,"_",diag,"_pheno=",pheno,"_MatrixEQTL_pheno=methylation_model=LINEAR.csv"),
        which_model="modelLINEAR"
      )
    }
  }
}

## run modelLINEAR_CROSS ##

# testing for the significance of the interaction between the genotype and the last covariate

for(cohort in c("PD", "PP", "BF")){
  
  run.matrixeqtl(
    SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_all_chrs_geno_matrix_012_wrangled-8b.csv"),
    expression_file_name = paste0("./MatrixEQTL_input/methylomics/",cohort,"_cohort_pheno=",pheno,"_methylation_matrix_MatrixEQTL_input.csv"),
    covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_cohort_cov_table_for_MatrixEQTL.csv"),
    output_file_name = paste0("./MatrixEQTL_output/",cohort,"_pheno=",pheno,"_cohort_MatrixEQTL_pheno=methylation_model=LINEARCROSS.csv"),
    which_model="modelLINEAR_CROSS"
  )
}
