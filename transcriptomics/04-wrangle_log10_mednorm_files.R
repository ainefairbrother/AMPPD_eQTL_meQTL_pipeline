# script that takes in log10 median normalised files (PP, PD, BF)
# filters for sample lists - which filters for visit month, gets rid of "Other" diagnoses, and filters for cohort
# prepares files for PEER

###### BF ###### 

bf.samples.to.use = read_csv("./data/AMPPD_sample_lists/BF-visit_month=0-case_control_other_latest!=Other_unique_1015.csv", show_col_types = FALSE) %>% 
  dplyr::pull(x)

bf = read.csv("./data/mt_nuc_corr_pipeline_output/BF_log10_mediannorm_TPM.csv", header=F) 
bf.new.colnames =  c("gene_id", bf[1,] %>% as.character()) %>% .[-length(.)]
colnames(bf) = bf.new.colnames

bf = bf %>% 
  tidyr::drop_na() %>% 
  tibble::column_to_rownames("gene_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate_if(is.character, as.numeric) %>% 
  `rownames<-`(gsub("\\.", "-", rownames(.))) %>% 
  .[rownames(.) %in% bf.samples.to.use, ] %>% 
  write_csv(x=., file="./data/PEER_input/BF_log10_mediannorm_TPM_wrangled.csv")

bf %>% dim()

###### PP ###### 

pp.samples.to.use = read_csv("./data/AMPPD_sample_lists/PP-visit_month=0-case_control_other_latest!=Other_unique_1015.csv", show_col_types = FALSE) %>% 
  dplyr::pull(x)

pp = read.csv("./data/mt_nuc_corr_pipeline_output/PP_log10_mediannorm_TPM.csv", header=F) 
pp.new.colnames =  c("gene_id", pp[1,] %>% as.character()) %>% .[-length(.)]
colnames(pp) = pp.new.colnames

pp = pp %>% 
  tidyr::drop_na() %>% 
  tibble::column_to_rownames("gene_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate_if(is.character, as.numeric) %>% 
  `rownames<-`(gsub("\\.", "-", rownames(.))) %>% 
  .[rownames(.) %in% pp.samples.to.use, ] %>% 
  write_csv(x=., file="./data/PEER_input/PP_log10_mediannorm_TPM_wrangled.csv")

###### PD ###### 

pd.samples.to.use = read_csv("./data/AMPPD_sample_lists/PD-visit_month=0-case_control_other_latest!=Other_unique_1015.csv", show_col_types = FALSE) %>% 
  dplyr::pull(x)

pd = read.csv("./data/mt_nuc_corr_pipeline_output/PD_log10_mediannorm_TPM.csv", header=F) 
pd.new.colnames =  c("gene_id", pd[1,] %>% as.character()) %>% .[-length(.)]
colnames(pd) = pd.new.colnames

pd = pd %>% 
  tidyr::drop_na() %>% 
  tibble::column_to_rownames("gene_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate_if(is.character, as.numeric) %>% 
  `rownames<-`(gsub("\\.", "-", rownames(.))) %>% 
  .[rownames(.) %in% pd.samples.to.use, ] %>% 
  write_csv(x=., file="./data/PEER_input/PD_log10_mediannorm_TPM_wrangled.csv")

###### all timepoint analysis ###### 

process.timepoint.files = function(samples.list.file.path, log.10.med.norm.file.path){
  
  out.file.name = gsub("TPM.csv", "TPM_wrangled.csv", log.10.med.norm.file.path) %>% 
    gsub("processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints", "PEER_input/all_timepoints", .)
  
  samples.to.use = read_csv(samples.list.file.path, show_col_types = FALSE) %>% 
    dplyr::pull(x) %>% 
    lapply(X=., FUN=function(X){return(str_match(string=X, pattern="P.+-")[1] %>% gsub("-$", "", .))}) %>% 
    unlist()
  
  log.10.med.norm = vroom::vroom(log.10.med.norm.file.path) 

  log.10.med.norm = log.10.med.norm %>% 
    tidyr::drop_na() %>% 
    tibble::column_to_rownames("gene_id") %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::mutate_if(is.character, as.numeric) %>% 
    `rownames<-`(gsub("\\.", "-", rownames(.))) %>% 
    .[rownames(.) %in% samples.to.use, ] %>% 
    write_csv(x=., file=out.file.name)
}

process.timepoint.files(samples.list.file.path="./data/AMPPD_sample_lists/PD-visit_month=ALL-case_control_other_latest!=Other_unique_1015.csv", 
                        log.10.med.norm.file.path="./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/PD_cohort_timepoint_averaged_log10_mediannorm_TPM.csv"
                        )
process.timepoint.files(samples.list.file.path="./data/AMPPD_sample_lists/PP-visit_month=ALL-case_control_other_latest!=Other_unique_1015.csv", 
                        log.10.med.norm.file.path="./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/PP_cohort_timepoint_averaged_log10_mediannorm_TPM.csv"
)
process.timepoint.files(samples.list.file.path="./data/AMPPD_sample_lists/PD-visit_month=ALL-case_control_other_latest==Case_unique_1015.csv", 
                        log.10.med.norm.file.path="./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/PD_case_timepoint_averaged_log10_mediannorm_TPM.csv"
)
process.timepoint.files(samples.list.file.path="./data/AMPPD_sample_lists/PD-visit_month=ALL-case_control_other_latest==Control_unique_1015.csv", 
                        log.10.med.norm.file.path="./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/PD_ctrl_timepoint_averaged_log10_mediannorm_TPM.csv"
)
process.timepoint.files(samples.list.file.path="./data/AMPPD_sample_lists/PP-visit_month=ALL-case_control_other_latest==Case_unique_1015.csv", 
                        log.10.med.norm.file.path="./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/PP_case_timepoint_averaged_log10_mediannorm_TPM.csv"
)
process.timepoint.files(samples.list.file.path="./data/AMPPD_sample_lists/PP-visit_month=ALL-case_control_other_latest==Control_unique_1015.csv", 
                        log.10.med.norm.file.path="./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/PP_ctrl_timepoint_averaged_log10_mediannorm_TPM.csv"
)







