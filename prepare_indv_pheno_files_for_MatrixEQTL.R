#' prepare.indv.pheno.files.for.MatrixEQTL
#'
#' @param file.path location to search for .csv files of the format rows=genes, cols=samples
#' @param file.pattern file pattern to select relevant files from the directory specified in file.path
#' @param loc.for.pheno.name.in.outfile this is the part of the input file name, filenamepart, that you wish to replace with filenamepart_pheno=
#'
#' @return returns .csv files for every row in the input file, labelled with the correct phenotype name 
#'
#' @examples

prepare.indv.pheno.files.for.MatrixEQTL = function(file.path, file.pattern, loc.for.pheno.name.in.outfile="log10"){
  
  print(paste(
    "Reading files from", file.path, "using pattern=", file.pattern
  ))
  
  file.list = list.files(path=file.path, pattern=file.pattern, full.names=FALSE)
  
  # define function to read in pheno files, split into single-row dfs and then write out into single-row df files
  read.file.and.split = function(x){
    
    f = read_csv(paste0(file.path,"/",x)) %>% 
      dplyr::filter(if_any(everything(), ~ !is.na(.)))
    
    f.split = f %>%
      tibble::rowid_to_column() %>%
      dplyr::group_split(rowid, .keep = FALSE) %>%
      setNames(., unique(f$sample_id))
    
    for(i in 1:length(f.split)){
      readr::write_csv(
        x=f.split[[i]], 
        file=paste0(file.path,"/",gsub(pattern=loc.for.pheno.name.in.outfile, 
                                       replacement=paste0("pheno=",names(f.split[i]), "_", loc.for.pheno.name.in.outfile),
                                       x=x)
        ))
    }
  }
  
  # apply function to file.list - all files in file.path with file.pattern
  lapply(X=file.list, FUN=read.file.and.split)
  
  print("File split done")
  
}

