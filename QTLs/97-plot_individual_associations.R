#' Title
#'
#' @param analysis_specific_dir_name if files are not within /MatrixEQTL_input/data-type/, and are within MatrixEQTL_input/transcriptomics/analysis_specific_dir_name instead, define this argument with the folder name
#' @param pheno ENS code for the MT gene from which expression is derived, or MT position for which the methylation prop. is derived
#' @param snp_chr_pos chr:pos format for SNP location
#' @param exp_or_meth looking at expression or methylation data
#' @param cohort PP, PD, BF (for single cohort) or PP_PD (for meta/mega)
#' @param diagnosis Case, Control (for LINEAR run), cohort (for LINEARCROSS run)
#'
#' @return
#' @export
#'
#' @examples
plot.indiv.assocs = function(analysis_specific_dir_name, pheno, snp_chr_pos, exp_or_meth, cohort, diagnosis, genomic_suffix=F){
  
  if(genomic_suffix==T){
    genomic_suffix = gsub(":","-", snp_chr_pos) %>% paste0("_chr",.)
  }
  
  snp_chr_pos = gsub(":", "-", snp_chr_pos) %>% paste0("chr", .)
  
  if(diagnosis=="cohort"){
    s1 = vroom::vroom(paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/",analysis_specific_dir_name,"/", cohort, "_all_chrs-maf0.05_geno_matrix_012_wrangled-8b",genomic_suffix,".csv"), show_col_types = FALSE) %>%
      tibble::column_to_rownames("pos") %>%
      .[snp_chr_pos, ]
    
  } else{
    s1 = vroom::vroom(paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/genomics/",analysis_specific_dir_name,"/", cohort, "_", tolower(diagnosis), "_all_chrs-maf0.05_geno_matrix_012_wrangled-8b",genomic_suffix,".csv"), show_col_types = FALSE) %>%
      tibble::column_to_rownames("pos") %>%
      .[snp_chr_pos, ]
  }
  if(exp_or_meth=="expression"){
    e1 = vroom::vroom(paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/transcriptomics/",analysis_specific_dir_name,"/", cohort, "_", diagnosis, "_pheno=",pheno,"_log10_mediannorm_TPM_masked_outliers_q1-3iqr_q3+3iqr_MT_ONLY.csv"), show_col_types = FALSE) %>%
      dplyr::filter(if_any(everything(), ~ !is.na(.))) %>%
      tibble::column_to_rownames("sample_id") %>%
      `colnames<-`(colnames(s1))
  }
  if(exp_or_meth=="methylation"){
    e1 = vroom::vroom(paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/methylomics/",analysis_specific_dir_name,"/", cohort, "_", diagnosis, "_pheno=", pheno,"_methylation_matrix_MatrixEQTL_input.csv"), show_col_types = FALSE) %>%
      dplyr::filter(if_any(everything(), ~ !is.na(.))) %>%
      tibble::column_to_rownames("sample_id") %>%
      `colnames<-`(colnames(s1))
    
    # e1[e1 == 1] = NA

  }
  
  cov = vroom::vroom(paste0("/home/abrowne/projects/amppd_analysis/data/MatrixEQTL_input/covariates/",analysis_specific_dir_name,"/", cohort, "_", diagnosis, "_maf0.05_cov_table_for_MatrixEQTL.csv"), show_col_types = FALSE) %>%
    tibble::column_to_rownames("variable") %>%
    `colnames<-`(colnames(s1))
  
  # bind all 3 
  all = rbind(e1, s1, cov) %>%
    t() %>%
    as.data.frame()
  
  formula = all[,pheno] ~ all[,snp_chr_pos]
  
  # + all[, "gPC1"] + all[, "gPC2"] + all[, "gPC3"] + all[, "gPC4"] + all[, "gPC5"] + all[, "PEER1"] + all[, "PEER2"] + all[, "PEER3"] + all[, "PEER4"] + all[, "PEER5"] + all[, "PEER6"] + all[, "PEER7"] + all[, "PEER8"] + all[, "PEER9"] + all[, "PEER10"] + all[, "sex"] + all[, "age_at_baseline"]
  
  lm1 = lm(formula=formula,
           data=all,
           na.action="na.exclude")
  
  # to allow vars to be used to select cols within ggplot call
  snp_chr_pos = sym(snp_chr_pos)
  pheno=sym(pheno)

  # generate counts of individuals in each genotype bin
  bin.n.counts = table(all[,2]) %>% 
    as.data.frame() 
  colnames(bin.n.counts) = c(snp_chr_pos, "Freq")

  print(bin.n.counts)
  
  return(
    all %>% 
      dplyr::filter(!!snp_chr_pos != -1) %>% # -1 denotes a missing value
      dplyr::filter(!!pheno != 1) %>%
      tidyr::drop_na() %>% 
      dplyr::mutate(geno.label=case_when(
        snp_chr_pos==0 ~ "AA",
        snp_chr_pos==1 ~ "AB",
        snp_chr_pos==2 ~ "BB"
      )) %>% 
      tibble::rownames_to_column("sample.id") %>% 
      dplyr::filter(sample.id %in% names(lm1$fitted.values)) %>% 

      ggplot(data=., aes(x=factor(!!snp_chr_pos), y=!!pheno)) +
      theme_minimal() +
      geom_boxplot() +
      geom_point(position="jitter", alpha=0.4) +
      geom_smooth(mapping=aes(group=1), method="lm", se=T, fullrange=F) +
      
      ggtitle(paste(cohort, diagnosis, pheno, snp_chr_pos)) +
      theme(plot.title = element_text(size = 8, face = "bold")) + 
      xlab("Genotype") +
      ylab(pheno)
  )
  
}