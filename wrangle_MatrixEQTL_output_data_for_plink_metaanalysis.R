
require(tidyverse)
require(vroom)
require(parallel)

# define function to help run another function in parallel across files in a dir
parallel.helper = function(fn, path, pattern, cores.to.use=3){
  
  # Run a function, fn, across files matching pattern, pattern, in directory, path. 
  # This will produce whatever files fn is designed to output
  file.list = list.files(path=path, pattern=pattern, full.names=T)
  parallel::mclapply(file.list, fn, mc.cores=cores.to.use)
  
}

wrangle.matrixeqtl.data.for.plink.metaanalysis = function(file.path){
  
  print(file.path)
  
  # if file size is less than 0.5MB, then this is an empty df, so don't process
  if( !(file.info(file.path)$size < 5e+05) ){
    
    vroom::vroom(file.path) %>% 
      # select relevant cols
      dplyr::select(SNP, beta, beta_se, p.value, snp.chr, snp.start) %>% 
      # rename to be compatible with plink --meta-analysis
      dplyr::rename(
        OR=beta,
        SE=beta_se,
        P=p.value,
        CHR=snp.chr,
        BP=snp.start
      ) %>%
      # add dummy cols for A1 and A2, as PLINK needs these cols to exist, but they're optional 
      # could map back to bed file to get them if necessary 
      dplyr::mutate(A1="A1",
                    A2="A2") %>% 
      # obtaining the exponent of a beta coefficient produces the OR - https://www.biostars.org/p/352366/
      dplyr::mutate(OR=exp(OR)) %>% 
      # write out as whitespace separated file to PLINK_input dir
      vroom::vroom_write(x=., path=gsub("MatrixEQTL_output", "PLINK_input", file.path) %>%
                           gsub("wrangled", "plinkIN", .), delim = " ")
  }
  
}

parallel.helper(fn=wrangle.matrixeqtl.data.for.plink.metaanalysis, 
                path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output", 
                pattern="wrangled.csv", 
                cores.to.use=4)
