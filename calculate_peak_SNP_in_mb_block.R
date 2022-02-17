assign.mb.block.to.group = function(group){
  # This function assigns, to a group of a df, block IDs, mb.block.id, to a cluster of row values
  # Rows are clustered based on the SNP position
  # Then the tree is cut to get positions that are within a series of megabase blocks
  # Positions in these blocks are given the same mb.block.id
  # It then determines the peak (min. p.value) SNP in the block 
  # The output is a modified group in tibble format 
  
  # define block size (bases)
  block.size=1000000
  
  # get distance tree by calculating the complete distance between positions
  tree = stats::hclust(stats::dist(group$snp.start), method = "complete")
  
  ##### test tier to check that the blocks (max-min) span 1mb  ##### 
  # split into 1MB clusters
  blocks = split(group$snp.start, stats::cutree(tree, h = block.size)) %>%
    `names<-`(paste0("block",names(.)))
  
  # check that the block range is less than 1mb - for each block
  error.count=0
  for(i in 1:length(blocks)){
    le.1mb = ((blocks[[i]] %>% max()) - (blocks[[i]] %>% min()) < block.size)
    if(le.1mb==FALSE){error.count=error.count+1}
  }
  if(error.count != 0){print(error.count)}
  
  #######################################################
  
  # split distance tree into 1MB clusters 
  return(
    split(group$snp.start, cutree(tree, h=block.size)) %>% 
      # name the clusters
      `names<-`(paste0("block",names(.))) %>% 
      # make list of 1MB blocks into df 
      stack() %>% 
      dplyr::rename(snp.start=values, mb.block.id=ind) %>% 
      # join to original df group
      dplyr::left_join(x=group, y=., by="snp.start")
  )
}
apply.assign.mb.fn.to.file.and.write.out = function(file.path){
  
  # This function takes in a file, and applies to it the assign.mb.block.to.group function which outputs a modified group
  # It does this using group_modify, which applies a function that outputs a modified group to all groups in a df
  # It then does some cleaning 
  # Then writes out the output, overwriting the input file 
  
  print(file.path)
  
  d = vroom::vroom(file.path) %>% 
    dplyr::select(-snp.chr) %>% 
    tidyr::extract(col=snp.chr.pos, into="snp.chr", regex="(chr.+):", remove=FALSE)
  
  # if the df is not empty, continue with wrangling
  if(dim(d)[1]!=0){

    # if the mb.block.id.chr is not present, continue with wrangling
    if( !("mb.block.id" %in% colnames(d)) ){
      d %>%
        dplyr::mutate(snp.chr=factor(snp.chr), 
                      phenotype=factor(phenotype), 
                      cohort=factor(cohort),
                      diagnosis=factor(diagnosis)) %>% 
        dplyr::group_by(snp.chr, phenotype, cohort, diagnosis) %>%
        # ensure group size is >=2, or there is nothing to cluster
        dplyr::filter(n() >= 2) %>% 
        dplyr::group_modify(~assign.mb.block.to.group(.x)) %>%
        dplyr::mutate(chr.mb.block.id=paste0(snp.chr,".",mb.block.id)) %>%
        # group by snp.chr and mb.block.id in order to get 'peak snp' in MB block - aka SNP with lowest p.value
        dplyr::group_by(chr.mb.block.id, phenotype, cohort, diagnosis, .drop=FALSE) %>%
        dplyr::mutate(min.pval.in.mb.block=min(p.value)) %>%
        dplyr::ungroup() %>%
        vroom::vroom_write(x=., file=file.path)
    }
  }
}
parallel.helper = function(fn, path, pattern, cores.to.use=3){
  
  # Run a function, fn, across files matching pattern, pattern, in directory, path. 
  # This will produce whatever files fn is designed to output
  file.list = list.files(path=path, pattern=pattern, full.names=T)
  parallel::mclapply(file.list, fn, mc.cores=cores.to.use)
  
}
lapply.helper = function(fn, path, pattern){
  
  # Run a function, fn, across files matching pattern, pattern, in directory, path.
  # This will produce whatever files fn is designed to output
  file.list = list.files(path=path, pattern=pattern, full.names=T)
  lapply(X=file.list, FUN=fn)
  
}
