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
    
    # # Assign rsid: prep data for convert_loc_to_rs.R
    # pos.to.snp = d %>% 
    #   dplyr::select(snp.chr, snp.start) %>% 
    #   dplyr::rename(CHR=snp.chr, BP=snp.start) %>% 
    #   # implement Regina's function from the coloc package - https://github.com/RHReynolds/colochelpR/blob/master/R/convert_rs_to_loc.R
    #   convert_loc_to_rs(df=., dbSNP = dbSNP) %>% 
    #   # rename to make cols compatible with matrixeqtl.out.df
    #   dplyr::rename(snp.chr=CHR, snp.start=BP) %>% 
    
    # # Assign rsid: re-join the original table to pos.to.snp - this is the table where pos are mapped to rsid 
    # dplyr::left_join(x=., 
    #                  y=pos.to.snp %>% 
    #                    tibble::tibble() %>% 
    #                    dplyr::mutate(snp.chr=as.character(snp.chr)) %>% 
    #                    dplyr::mutate(snp.start=as.numeric(snp.start)), 
    #                  by=c("snp.chr", "snp.start")) %>% 
    # dplyr::distinct()
  }
  parallel::mclapply(file.list, wrangle.write.file, mc.cores=4)
}