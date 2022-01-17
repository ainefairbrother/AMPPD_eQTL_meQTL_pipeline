#!/usr/bin/env Rscript

# import libs
require(data.table)
require(tidyverse)

wrangle.geno.mat.for.MatrixEQTL = function(pos.file, geno.matrix.file, sample.id.file, out.file){
  
  # import pos and make character col label vector by pasting chr and numeric position together
  # import the pos (genomic pos) file
  print(paste("Importing", pos.file))
  pos=data.table::fread(pos.file, header=F) %>% 
    dplyr::mutate(V3 = paste0(V1, "-", V2)) %>% 
    dplyr::pull(V3)

  # import mat
  print(paste("Importing", geno.matrix.file))
  mat = data.table::fread(geno.matrix.file)
  
  # import sample id list - samples that have transcriptomics, methylomics, genomics & appropriate metadata
  samples.to.use = readr::read_csv(sample.id.file) %>% 
    dplyr::pull(".")
  
  if(cohort=="BF"){samples.to.use = gsub("-SVM0_5T1", "", samples.to.use)}
  if( (cohort=="PP") | (cohort=="PD") ){samples.to.use = gsub("-BLM0T1", "", samples.to.use)}
  
  # check dim
  mat %>% dim() %>% print()
  
  # check head
  mat[1:5,1:5]
  
  # check mat and pos dims
  pos %>% length() %>% print()
  mat %>% dim() %>% print()
  
  # filter columns for samples.to.use
  mat = mat[ ,colnames(mat) %in% samples.to.use, with=FALSE]
  
  # assign pos as rownames to the geno matrix
  mat$pos = pos
  
  # check mat
  print("Checking rownames, colnames, dim and head")
  rownames(mat)[1:5] %>% print()
  colnames(mat)[1:5] %>% print()
  dim(mat) %>% print()
  mat[1:10,1:10] %>% print()
  
  # filter and write out
  print("Writing data...")
  
  fwrite(x=mat %>% dplyr::relocate(pos), file=out.file)
  
  rm(mat)
  
}

for(cohort in c("PD", "PP", "BF")){
base.dir="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/"
  for(diag in c("case", "control")){

    if(diag=="case"){diag.sample="Case"}
    if(diag=="control"){diag.sample="Control"}

    wrangle.geno.mat.for.MatrixEQTL(
      pos.file = paste0(base.dir,cohort,"_all_chrs-maf0.05_geno_matrix.012.pos"),
      geno.matrix.file = paste0(base.dir,cohort,"_",diag.sample,"_all_chrs-maf0.05_geno_matrix_012_wrangled.csv"),
      sample.id.file = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/",cohort,"_",diag.sample,"_reduce_intersect_gpc_peer_meta_sampleid.csv"),
      out.file = paste0(base.dir,cohort,"_",diag,"_all_chrs-maf0.05_geno_matrix_012_wrangled-8b.csv")
    )

  }
}

for(cohort in c("PD", "PP", "BF")){

  base.dir="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/"

  wrangle.geno.mat.for.MatrixEQTL(
    pos.file = paste0(base.dir,cohort,"_all_chrs-maf0.05_geno_matrix.012.pos"),
    geno.matrix.file = paste0(base.dir,cohort,"_cohort_all_chrs-maf0.05_geno_matrix_012_wrangled.csv"),
    sample.id.file = paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/",cohort,"_cohort_reduce_intersect_gpc_peer_meta_sampleid.csv"),
    out.file = paste0(base.dir,cohort,"_all_chrs-0.05_geno_matrix_012_wrangled-8b.csv")
  )
}

