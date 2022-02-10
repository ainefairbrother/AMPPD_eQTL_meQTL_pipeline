# this script is to filter the PP-celltype sample covariate files, to remove the PEER factor rows,
# producing covariate files with the standard covs + celltype covs 
# this is to check whether colinearity of the PEER and celltype covariates is causing the 
# in/de -flation of p.values observed

in.dir = "/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/celltype_plus_standard_correction_with_PEER//"
out.dir = "/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/celltype_plus_standard_correction_no_PEER/"

# get files in the celltype_PP_samples_filtered dir 
files.to.edit = list.files(path=in.dir, pattern=".csv")

# edit the files to remove PEER factor rows and save in the celltype_plus_standard_correction_no_PEER dir
for(f in files.to.edit){
  print(f)
  
  readr::read_csv(file=paste0(in.dir,f)) %>% 
    dplyr::filter(!grepl(pattern="PEER", x=variable)) %>% 
    write_csv(x=., file=paste0(out.dir, f))

}

