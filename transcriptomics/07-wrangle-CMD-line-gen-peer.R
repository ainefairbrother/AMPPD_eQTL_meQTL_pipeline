# this script assigns sample IDs to the PEER factor output files which are returned from the CMD line PEER tool without labels
# load the PEER input file, grab the row names, and assign them to the PEER output file 

for(fname in c("PP_case", "PP_ctrl", "PD_case", "PD_ctrl", "PP_cohort", "PD_cohort")){
  
  rn = vroom(paste0("./data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints/",fname,"_timepoint_averaged_log10_mediannorm_TPM.csv")) %>% 
    dplyr::select(-gene_id) %>% 
    colnames(.)
  
  dat = read.csv(paste0("./data/PEER_output/rosalind-cmd-line-gen/",fname,"_timepoint_averaged_PEER.csv"), header = F) %>% 
    t() %>% 
    `rownames<-`(c(rn)) %>% 
    `colnames<-`(c(paste0("V", c(1:50))))
  
  dat %>% 
    write.csv(x=., paste0("./data/PEER_output/rosalind-cmd-line-gen/",fname,"_timepoint_averaged_PEER_wrangled.csv"))
  
}