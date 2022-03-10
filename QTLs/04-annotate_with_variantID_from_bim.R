#' Title
#'
#' @param target.table 
#' @param refalt 
#'
#' @return
#' @export
#'
#' @examples
annotate.CHR.BP.table.with.variantid = function(target.table, ref.alt=T, diff.cohort=NULL){
  
  library(genio)
  
  # read in bim files
  bim = dplyr::bind_rows(
    genio::read_bim(file="/home/abrowne/projects/amppd_analysis/data/genomics_bed_bim_fam_from_rosalind/PD_all_chrs_maf0.05.bim", verbose = TRUE) %>% dplyr::mutate(cohort="PD"),
    genio::read_bim(file="/home/abrowne/projects/amppd_analysis/data/genomics_bed_bim_fam_from_rosalind/PP_all_chrs_maf0.05.bim", verbose = TRUE) %>% dplyr::mutate(cohort="PP")
    )
  
  # if diff.cohort arg is supplied as a character, as opposed to NULL, use this to generate the label column 'cohort' instead of matching the cohort in the file name
  # this is for joint analyses like mega or meta, where the associations are derived from combined cohort samples 
  # primarily allows this function to map variant IDs to CHR-BP in both tables
  if(!is.null(diff.cohort)){
    bim=bim %>% 
      dplyr::mutate(cohort=diff.cohort) %>% 
      dplyr::distinct()
  }
  
  if(ref.alt==T){
    bim = bim %>% 
      # make variant.id column
      dplyr::mutate(variant.id=paste0(chr, "_", pos, "_", ref, "_", alt))
  }
  if(ref.alt==F){
    bim = bim %>% 
      # make variant.id column
      dplyr::mutate(variant.id=paste0(chr, "_", pos, "_", alt, "_", ref))
  }
  
  # add relevant variant IDs to target.table
  target.table %>%
    # ensure joining col exists
    dplyr::mutate(snp.chr.pos=paste0("chr",snp.chr,":",snp.start)) %>% 
    dplyr::left_join(x=., y=bim %>%
                       dplyr::mutate(snp.chr.pos=paste0("chr",chr,":",pos)) %>%
                       dplyr::select(snp.chr.pos, variant.id, cohort),
                     by=c("snp.chr.pos", "cohort")) %>%
    return(.)
}