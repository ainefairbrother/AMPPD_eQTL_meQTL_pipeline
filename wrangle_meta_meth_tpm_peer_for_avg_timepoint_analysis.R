#!/usr/bin/env Rscript

# run with: 
# screen -S run-avg-tpm
# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/
# Rscript wrangle_meta_meth_tpm_peer_for_avg_timepoint_analysis.R

##################### load libraries ######################
require(dplyr)
require(tidyr)
require(vroom)

##################### apply to cohort-diagnosis files ######################

calc.avg.tpm.across.timepoints = function(f){
  
  print(paste("Processing ", f))
  
  vroom::vroom(f, show_col_types = FALSE) %>% dim() %>% print()
  
  vroom::vroom(f, show_col_types = FALSE) %>% 
    tidyr::pivot_longer(cols=2:ncol(.), names_to="sample_id", values_to="log10_mednorm_TPM") %>% 
    dplyr::rename("gene_id"=`...1`) %>% 
    dplyr::mutate(sample_id=gsub("\\.", "-", sample_id)) %>% 
    tidyr::extract(col=sample_id, into="participant_id", regex="(P.+)-", remove=F) %>% 
    
    # get timepoint mean of TPMs
    dplyr::group_by(participant_id, gene_id) %>% 
    dplyr::mutate(mean.log10.mednorm.tpm = mean(log10_mednorm_TPM, na.rm=T)) %>% 
    
    # ungroup and remove cols that are no longer needed
    dplyr::ungroup() %>% 
    dplyr::select(-c(log10_mednorm_TPM, sample_id)) %>% 
    
    # need to make cols unique, because now there is a repeat row for every timepoint that previously existed
    dplyr::distinct() %>% 
    
    # pivot back into matrix form 
    tidyr::pivot_wider(id_cols=c("gene_id"),
                       values_from="mean.log10.mednorm.tpm", 
                       names_from="participant_id") %>% 
    
    # write out
    vroom::vroom_write(x=., 
                       file=paste0(gsub("_all_timepoint_log10_mediannorm_TPM.csv", "_timepoint_averaged_log10_mediannorm_TPM.csv", f)))
  
  vroom::vroom(paste0(gsub("_all_timepoint_log10_mediannorm_TPM.csv", "_timepoint_averaged_log10_mediannorm_TPM.csv", f)), show_col_types = FALSE) %>% dim() %>% print()
  
  print("done")
  
}

file.list = list.files(path="/home/abrowne/projects/amppd_analysis/data/processed_and_unprocessed_transcriptomic/mt_nuc_corr_pipeline_output/all_timepoints", 
                       pattern="_all_timepoint_log10_mediannorm_TPM.csv", full.names = T)

lapply(X=file.list, FUN=calc.avg.tpm.across.timepoints)

