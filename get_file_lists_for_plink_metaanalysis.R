get.file.lists.for.plink.metaanalysis = function(file.path, file.suffix="plinkIN.csv", cohort1="PP", cohort2="PD"){
  
  cohort1.pattern = paste0(cohort1, ".+", file.suffix)
  cohort2.pattern = paste0(cohort2, ".+", file.suffix)
  
  cohort1.files = list.files(path=file.path, pattern=cohort1.pattern, full.names=F)
  cohort2.files = list.files(path=file.path, pattern=cohort2.pattern, full.names=F)
  
  # get diagnosis-pheno-model files in common between the two datasets
  # remove cohort from file name and get intersecting diagnosis-pheno-model files
  intersect(
    cohort1.files %>% gsub("^[A-Z]{2}", "", .), 
    cohort2.files %>% gsub("^[A-Z]{2}", "", .) 
  ) %>% 
    # write to plain text file
    writeLines(., file(paste0(file.path, "/", cohort1, "_", cohort2, "_file_list_metaanalysis_plinkIN.txt")), sep="\n")
}

get.file.lists.for.plink.metaanalysis(
  file.path="/home/abrowne/projects/amppd_analysis/data/PLINK_input",
  file.suffix="plinkIN.csv",
  cohort1="PP",
  cohort2="PD"
)