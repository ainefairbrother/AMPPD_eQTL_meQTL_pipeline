run.matrixeqtl = function(
  SNP_file_name, 
  pheno_file_name, 
  covariates_file_name,
  which_pheno, # can be expression or methylation 
  output_file_name,
  which_model
){
  
  # # test
  # cohort = "PP"
  # diag = "Control"
  # if(diag=="Case"){diag.lower="case"}
  # if(diag=="Control"){diag.lower="control"}
  # 
  # SNP_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/all_timepoints/test.csv")
  # # include pheno= in pheno_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  # pheno_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/transcriptomics/all_timepoints/",cohort,"_",diag,"_timepoint_averaged_pheno=_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv")
  # covariates_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/all_timepoints/",cohort,"_",diag,"_averaged_timepoint_maf0.05_cov_table_for_MatrixEQTL.csv")
  # # include pheno= in output_file_name as a placemarker for the run.matrixeqtl function to replace it with the actual phenotype label
  # output_file_name = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/all_timepoints/",cohort,"_",diag,"_pheno=_maf=0.05_phenotype=expression_model=LINEAR_MatrixEQTL.csv")
  # which_model="modelLINEAR"
  # which_pheno="expression"
  
  # load packages
  # MatrixEQTL package - http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
  library(MatrixEQTL)
  
  # set working dir
  setwd("/home/abrowne/projects/amppd_analysis/data/")
  
  # Then set the parameters such as selected linear model and names of genotype and expression data files.
  # modelANOVA or modelLINEAR or modelLINEAR_CROSS, here the genotype is assumed to have only additive effect on expression
  # additive will test for linear additive model - i.e you will compare AA vs AB vs BB in an additive way (so the SNP allele B has twice the effect when present in two copies)
  # the dominant model will group inds with at least one alternative allele, so AA vs (AB or BB)
  # if you think effects are dominant you could use this approach, or if you have too few samples and you want to increase power
  # the additive model is probably best, as we have a reasonable sample number to jusitfy 3 bins 
  if(which_model=="modelLINEAR"){useModel=modelLINEAR}
  if(which_model=="modelANOVA"){useModel=modelANOVA}
  if(which_model=="modelLINEAR_CROSS"){useModel=modelLINEAR_CROSS}
  
  # define locations of genotype matrix, expression/methylation matrix, covariate file, output file
  # A separate file may be provided with extra covariates. 
  # In case of no covariates set the variable covariates_file_name to character().
  # file names included in function args:
  # SNP_file_name
  # pheno_file_name
  # covariates_file_name
  # output_file_name
  
  # The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. 
  # Note that for larger datasets the threshold should be lower. 
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis > 0 to perform eQTL analysis with separate p-value 
  # thresholds for local and distant eQTLs. Distant and local associations significant at corresponding 
  # thresholds are recorded in output_file_name and output_file_name.cis respectively and in the returned object. 
  # In this case the false discovery rate is calculated separately for these two sets of eQTLs.
  pvOutputThreshold = 0.05 # simple cut-off to minimise output-file size
  
  # Finally, define the covariance matrix for the error term.
  # This parameter is rarely used. 
  # If the covariance matrix is a multiple of identity, set it to numeric().
  errorCovariance = numeric()
  
  # The next section of the sample code contains three very similar parts loading the files with 
  # genotype, gene expression, and covariates. 
  # In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), 
  # the string representation for missing values, the number of rows with column labels, 
  # and the number of columns with row labels. Finally, one can change the number of the variables
  # in a slice for the file reading procedure (do not change if not sure).
  
  # load genotype matrix
  snps = SlicedData$new()
  snps$fileDelimiter = ","       # define the sep of the file - i.e comma or space separated file
  snps$fileOmitCharacters = "-1"  # denote missing values - the missing values are defined as -1 by vcftools --012
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # one column of row labelsjavascript:
  snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  snps$LoadFile(SNP_file_name)
  
  # load covariate matrix 
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = ","       # comma separated file
  cvrt$fileOmitCharacters = "NA" # denote missing values
  cvrt$fileSkipRows = 1          # one row of column labels
  cvrt$fileSkipColumns = 1       # one column of row labels
  cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  cvrt$LoadFile(covariates_file_name)
  
  # define phenos, an object containing a vector of all the individual phenotypes in the phenotype file, 
  # this will be genes for the expression data, and MT positions for the methylation data
  if(which_pheno=="expression"){phenos=c("ENSG00000198888", "ENSG00000198763", "ENSG00000198840",
                                         "ENSG00000212907", "ENSG00000198886", "ENSG00000198786", 
                                         "ENSG00000198695","ENSG00000198712", "ENSG00000198938", 
                                         "ENSG00000228253", "ENSG00000198727", "ENSG00000198804", 
                                         "ENSG00000198899", "ENSG00000211459","ENSG00000210082")}
  
  if(which_pheno=="methylation"){phenos=c("3238_full", "4271_full", "4392_full", "5520_full", 
                                          "5647_full", "5721_full", "5818_full", "5883_full", 
                                          "7526_full", "8303_full", "9999_full", "10413_full", 
                                          "12146_full", "12274_full", "14734_full", "15896_full", 
                                          "15948_full", "2617_full", "13710_full")}
  
  # loop through the individual phenotypes, read in pheno file (into obj called gene), run matrix eqtl and produce a qq plot for QC purposes
  for(pheno in phenos){
    
    # generate pheno input file name and output file name
    indv.pheno.file.name = gsub("pheno=", paste0("pheno=", pheno), pheno_file_name)
    indv.output.file.name = gsub("pheno=", paste0("pheno=", pheno), output_file_name)
    
    # test whether pheno file exists, if it does not, skip to the next phenotype 
    if(file.exists(indv.pheno.file.name)){
      
      print(paste("Generating", indv.output.file.name))
      
      # load expression matrix
      gene = SlicedData$new()
      gene$fileDelimiter = ","       # comma separated file
      gene$fileOmitCharacters = "NA" # denote missing values
      gene$fileSkipRows = 1          # one row of column labels
      gene$fileSkipColumns = 1       # one column of row labels
      gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
      gene$LoadFile(indv.pheno.file.name)
      
      # Finally, the main Matrix eQTL function is called:
      me = Matrix_eQTL_engine(
        snps = snps, # SNPs
        gene = gene, # pheno
        cvrt = cvrt, # covariates
        output_file_name = indv.output.file.name, # name of the output file
        # Set pvOutputThreshold > 0 and pvOutputThreshold.cis = 0 (or use Matrix_eQTL_engine)
        # to perform eQTL analysis without using gene/SNP locations. Associations significant at the
        # pvOutputThreshold level are be recorded in output_file_name and in the returned object
        pvOutputThreshold = pvOutputThreshold, # Significance threshold for all/distant tests.
        useModel = useModel, # which model to use, define above
        errorCovariance = errorCovariance,
        verbose = FALSE, # get progress report
        min.pv.by.genesnp = FALSE, # record the minimum p-value for each SNP and each gene in the returned object                                                                                                                                                        
        noFDRsaveMemory = FALSE,
        pvalue.hist = "qqplot")
      
      # generate qq plot and save
      tiff(gsub(".csv", "_qq.tiff", indv.output.file.name)) 
      plot(me, pch = 16, cex = 0.7)
      dev.off() 
      
      # Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me. 
      # In case of cis/trans eQTL analysis described below, two output files are produced, one with cis-eQTLs, another only with trans. 
      # Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.
    } else{next}
  }
}


