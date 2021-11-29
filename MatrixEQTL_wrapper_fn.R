
run.matrixeqtl = function(
  SNP_file_name, 
  expression_file_name, 
  covariates_file_name,
  output_file_name,
  which_model="modelLINEAR"
){
  
  # test
  # cohort = "BF"
  # diag = "Case"
  # diag.lower = "case"
  # SNP_file_name = paste0("./MatrixEQTL_input/genomics/",cohort,"_",diag.lower,"_all_chrs_geno_matrix_012_wrangled-8b.csv")
  # expression_file_name = paste0("./MatrixEQTL_input/transcriptomics/",cohort,"_",diag,"_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv")
  # covariates_file_name = paste0("./MatrixEQTL_input/covariates/",cohort,"_",diag,"_cov_table_for_MatrixEQTL.csv")
  # output_file_name = paste0("./MatrixEQTL_output/",cohort,"_",diag,"_MatrixEQTL_pheno=expression_model=LINEAR.csv")
  # which_model="modelLINEAR"
  
  # load the package - http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
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
  
  # define locations of genotype matrix and expression matrix
  SNP_file_name = SNP_file_name
  expression_file_name = expression_file_name
  
  exp.file = read_csv(expression_file_name)

  
  # A separate file may be provided with extra covariates. In case of no covariates set the variable covariates_file_name to character().
  covariates_file_name = covariates_file_name
  
  # define the outfile loc and name
  output_file_name = output_file_name
  
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
  snps$fileDelimiter = ","       # comma separated file
  snps$fileOmitCharacters = "NA" # denote missing values
  snps$fileSkipRows = 1          # one row of column labels
  snps$fileSkipColumns = 1       # one column of row labelsjavascript:
  snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  snps$LoadFile(SNP_file_name)
  
  # load expression matrix
  gene = SlicedData$new()
  gene$fileDelimiter = ","       # comma separated file
  gene$fileOmitCharacters = "NA" # denote missing values
  gene$fileSkipRows = 1          # one row of column labels
  gene$fileSkipColumns = 1       # one column of row labels
  gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  #gene$LoadFile(expression_file_name)

  # load covariate matrix 
  cvrt = SlicedData$new()
  cvrt$fileDelimiter = ","       # comma separated file
  cvrt$fileOmitCharacters = "NA" # denote missing values
  cvrt$fileSkipRows = 1          # one row of column labels
  cvrt$fileSkipColumns = 1       # one column of row labels
  cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
  cvrt$LoadFile(covariates_file_name)
  
  # Finally, the main Matrix eQTL function is called:
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    pvalue.hist = "qqplot")
  
  # generate qq plot and save
  tiff(gsub(".csv", "_qq.tiff", output_file_name)) 
  plot(me, pch = 16, cex = 0.7)
  dev.off() 
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    pvalue.hist = 100)
  
  # generate pval dist plot and save
  tiff(gsub(".csv", "_pdist.tiff", output_file_name)) 
  plot(me, col="grey")
  dev.off() 
  
  # save .me object 
  save(me, file=gsub(".csv", ".Rdata", output_file_name))
  
  # Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me. 
  # In case of cis/trans eQTL analysis described below, two output files are produced, one with cis-eQTLs, another only with trans. 
  # Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.
  
}


