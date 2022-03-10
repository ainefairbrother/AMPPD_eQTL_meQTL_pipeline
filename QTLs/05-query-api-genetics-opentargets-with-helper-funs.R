#### FUNS #####

#' query.opentargets.with.varIDs
#'
#' @param variantID This function queries the open targets genetics API with a variant ID in the format: chr_pos_alt_ref i.e. 1_154453788_C_T
#'
#' @return This function returns a df with 1 row corresponding to the query entered. 
#' The output can then be applied to a vector containing variant_IDs and then the rows bound to generate a df containing the following cols: 
#' variant_id, chromosome, position, refAllele, altAllele, nearestGene.id, nearestGene.symbol, nearestGene.bioType
query.opentargets.with.varIDs = function(variantId){
  
  # Import libs
  library(httr)
  library(tidyverse)
  
  # Build query string
  query_string="
  query var($variantId: String!){
  variantInfo(variantId: $variantId){
  id
  chromosome
  position
  refAllele
  altAllele
  nearestGene {
  id
  symbol
  bioType
  }
  }
  }"
  
  # Set base URL of Genetics Portal GraphQL API endpoint
  base_url <- "https://api.genetics.opentargets.org/graphql"
  
  # Set variables object of arguments to be passed to endpoint
  variables <- list("variantId" = variantId)
  
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)
  
  # Perform POST request
  r <- httr::POST(url=base_url, body=post_body, encode='json')
  
  # make res (list of lists) into df
  r.data.slot = content(r)$data 
  
  if( is.null(r.data.slot$variantInfo$nearestGene)==TRUE ){
    r.data.slot$variantInfo$nearestGene = NA
  }
  
  return(r.data.slot %>%
           lapply(., data.frame, stringsAsFactors = FALSE) %>%
           dplyr::bind_rows(.))
}

# define function to apply query.opentargets.with.varIDs to a list of variantIDs
#' apply.to.variantID.list.and.generate.output.table
#' This is a helper function to run the `query.opentargets.with.varIDs` function on a list of variantIDs in parallel and produce a single output dataframe. 
#'
#' @param varID.list This is a list of variantIDs in the format: chr_pos_alt_ref i.e. 1_154453788_C_T. 
#' @param cores.to.use Number of cores to use to run the query - parellisation helps to speed this up, as inetraction with the API for large lists is slow. 
#'
#' @return This function returns a df with nrow=length(varID.list) corresponding to the query entered. 
#' The output can then be applied to a vector containing variant_IDs and then the rows bound to generate a df containing the following cols: 
#' variant_id, chromosome, position, refAllele, altAllele, nearestGene.id, nearestGene.symbol, nearestGene.bioType
#' @example `query.opentargets.with.varIDs.helper(varID.list=c("1_154453788_C_T", "1_154453788_C_A"), cores.to.use=2)`  
#' 
query.opentargets.with.varIDs.helper = function(varID.list, cores.to.use=2){
  
  library(tidyverse)
  library(parallel)
  
  mclapply(X=varID.list, FUN=query.opentargets.with.varIDs, mc.cores=cores.to.use) %>%
    dplyr::bind_rows() %>% 
    tibble::as_tibble() %>% 
    return()
  
}

# define function to apply query.opentargets.with.varIDs to all variantIDs in a file
apply.opentarget.query.to.variantid.column.FILE = function(file.path, mc.cores=10){
  
  library(parallel)
  library(dplyr)
  library(tibble)
  library(vroom)
  
  in.df = vroom::vroom(file=file.path) %>% 
    data.table::as.data.table()
  
  if("variant.id" %in% colnames(in.df)){
    parallel::mclapply(X=in.df$variant.id, FUN=query.opentargets.with.varIDs, mc.cores=mc.cores) %>%
      dplyr::bind_rows() %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(snp.chr.pos=paste0("chr",chromosome,":",position)) %>%
      dplyr::select(any_of( c("snp.chr.pos", "nearestGene.id", "nearestGene.symbol", "nearestGene.bioType") )) %>%
      dplyr::distinct() %>%
      dplyr::left_join(x=in.df , y=., by="snp.chr.pos") %>% 
      vroom::vroom_write(x=., file=file.path)
  }
}

# define function to apply query.opentargets.with.varIDs to all variantIDs in a df
apply.opentarget.query.to.variantid.column.DF = function(target.table, mc.cores=10){
  
  library(parallel)
  library(dplyr)
  library(tibble)
  
  if("variant.id" %in% colnames(target.table)){
    parallel::mclapply(X=target.table$variant.id, FUN=query.opentargets.with.varIDs, mc.cores=mc.cores) %>%
      dplyr::bind_rows() %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(snp.chr.pos=paste0("chr",chromosome,":",position)) %>%
      dplyr::select(any_of( c("snp.chr.pos", "nearestGene.id", "nearestGene.symbol", "nearestGene.bioType") )) %>%
      dplyr::distinct() %>%
      dplyr::left_join(x=target.table , y=., by="snp.chr.pos") %>% 
      return(.)
  }
}

loop.helper = function(fn, path, pattern){
  
  file.list = list.files(path=path, pattern=pattern, full.names=T) 
  for(i in 1:length(file.list)){
    fn(file.list[i])
  }
}

# #### Run on files in path #####
# 
# loop.helper(fn=apply.to.df.and.wrangle, 
#             path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/", 
#             pattern="P.+_wrangled.csv")

