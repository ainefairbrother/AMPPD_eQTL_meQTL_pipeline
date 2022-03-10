#!/usr/bin/env Rscript

# Author: Aine Fairbrother-Browne
# Date: 18/22
# Description: 
# cd /home/abrowne/projects/amppd_analysis/pipeline/eqtls_mqtls
# conda activate r4-base
# Rscript 99-join_cohort_mega_meta.R

#### -------------- prepare env -------------- ####

# conda activate r4-base
# screen -S 
# conda activate r4-base

#### -------------- load libs -------------- ####
require(dplyr)
require(tidyr)
require(vroom)

#### -------------- define functions -------------- ####

prepare.res.table.for.joint = function(res.table, unique.col.suffix){
  res.table %>% 
    dplyr::filter(cohort!="BF") %>% # the individual cohort has BF, so remove
    dplyr::select_all(~gsub("p.value", "P", .)) %>% # rename p.value to P to match mega analysis table which is a PLINK output - this is the fixed effects P value
    dplyr::select(-any_of(c("nuc.snp.maf.filter", "phenotype.category", "t.stat", "beta_se", "snp.chr.pos", "beta", "FDR"))) %>% # meta doesn't have a beta, it has an OR, so remove this col until I find a way to harmonise
    dplyr::select(-any_of(c("P.R."))) %>% # this is the P value for random effects in the meta analysis, so remove this for the meta table
    dplyr::select(-contains("block"), -contains("nearestGene")) %>% 
    dplyr::rename_at(vars(P), function(x) paste0(x,unique.col.suffix)) %>% # add suffix to col names that I want to keep to make them unique to the df being passed in - will be needed when the tables are joined together
    dplyr::select_all(~gsub("snp.start", "BP", .)) %>% # rename snp.start as BP
    dplyr::select_all(~gsub("snp.chr", "CHR", .)) %>% # rename snp.chr as CHR, to harmonise with PLINK mega-analysis output, the select_all(~gsub()) method is useful because it is conditional i.e. only does it if the col contains snp.chr
    #dplyr::mutate(CHR=gsub("chr","",CHR)) %>% # remove chr from chromosome labels in CHR col
    dplyr::mutate(unique.assoc.id = paste0(cohort,"-",diagnosis,"-",matrixeqtl.model.run,"-",CHR,"-",BP)) %>% # make unique id
    dplyr::select(unique.assoc.id, contains("P_")) %>% # 
    #dplyr::relocate(unique.assoc.id, phenotype, BP, contains("P")) %>% 
    return(.)
}
join.analyses.write.out = function(cohort.path, mega.path, meta.path, out.path){
  
  dplyr::full_join(
    # apply prepare.res.table.for.joint to harmonise and then join meta and mega by unique.assoc.id
    x=vroom::vroom(mega.path, show_col_types = FALSE) %>% 
      tidyr::extract(col=snp.chr.pos, into="snp.chr", regex="chr(.+):", remove=FALSE) %>% 
      prepare.res.table.for.joint(res.table=., unique.col.suffix="_mega"),
    y=vroom::vroom(meta.path, show_col_types = FALSE) %>% 
      prepare.res.table.for.joint(res.table=., unique.col.suffix="_meta"),
    by=c("unique.assoc.id")
  ) %>% 
    dplyr::arrange(P_mega, P_meta) %>% 
    # join the mega-meta columns to the cohort table
    dplyr::full_join(x=., y=vroom::vroom(cohort.path, show_col_types = FALSE) %>% 
                       # tidyr::extract(col=snp.chr.pos, into="snp.chr", regex="chr(.+):", remove=FALSE) %>% 
                       dplyr::rename(CHR=snp.chr, BP=snp.start) %>% 
                       dplyr::mutate(unique.assoc.id = paste0("PP_PD","-",diagnosis,"-",matrixeqtl.model.run,"-",CHR,"-",BP)), # add unique id to cohort res, adding PP_PD to the start to match up with mega and meta, review this??
                     by="unique.assoc.id") %>% 
    dplyr::filter(!is.na(p.value)) %>% 
    vroom::vroom_write(x=., file=out.path)
  
}

#### -------------- implement functions -------------- ####

join.analyses.write.out(
  cohort.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/aggregated_tables/matrixeqtl_res_aggregated_no_p_filter.csv", 
  mega.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/aggregated_tables/matrixeqtl_mega_res_aggregated_no_p_filter.csv", 
  meta.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/aggregated_tables/matrixeqtl_meta_res_aggregated_no_p_filter.csv", 
  out.path="/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_output/aggregated_tables/matrixeqtl_res_aggregated_no_p_filter_cohort_mega_meta.csv"
)