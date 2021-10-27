
## 1. Setup

library(bigsnpr)
library(tidyverse)

setwd("/home/abrowne/projects/amppd_analysis/")
plot.out="/home/abrowne/projects/amppd_analysis/results/AMPPD_gPC_results/"
data.out="/home/abrowne/projects/amppd_analysis/data/bigsnpr_gPC_output/"

# set up analysis loop
for(cohort in c("BF", "PP", "PD")){
  
  print("----------------------")
  print(paste("Running analysis for cohort=", cohort))
  
  ## 2. Load data
  
  # Get the example bedfile from package bigsnpr
  bedfile <- paste0("./data/genomics_bed_bim_fam_from_rosalind/",cohort,"_all_chrs.bed")
  
  ## 3. Pre-prepare data
  
  # Read from bed/bim/fam, it will create new files - let's put them in a temporary directory (tmpfile) for this demo.
  tmpfile <- tempfile()
  
  snp_readBed(bedfile, backingfile = tmpfile)
  
  # Attach the "bigSNP" object in R session
  obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))
  
  # See what it looks like
  # str(obj.bigSNP, max.level = 2, strict.width = "cut")
  
  # Get aliases for useful obj.bigSNP slots
  CHR <- obj.bigSNP$map$chromosome %>% 
    # getting the following error: "Error: 'infos.chr' should contain only integers.", so need to convert X and Y into numerics
    # so assign chrX as chr23 and chrY as chr24
    replace(., .=="X", 23) %>%
    replace(., .=="Y", 24) %>%
    as.numeric()
  
  POS <- obj.bigSNP$map$physical.pos
  
  print("Running imputation")
  
  G <- obj.bigSNP$genotypes %>% 
    # svd will not accept missing values in the genotype matrix so running imputation
    # snp_fastImpute function benchmarked against widely used imputation software beagle 
    # found that beagle made 3.1% of imputation errors, snp_fastImpute made 4.7%
    # but snp_fastImpute was 20X faster
    # snp_fastImpute made less 0/2 switching errors (imputing with a homozygous ref where the true genotype is homozygous variant)
    # snp_fastImpute was throwing xgboost errors, so utilising snp_fastImputeSimple instead
    snp_fastImputeSimple(Gna=., method = "random") # random = sampling according to allele frequencies
   
  # Check some counts for the 10 first SNPs
  big_counts(G, ind.col = 1:50)
  
  print("Generating svd0 object")
  
  ## 4. Run genetic PCA
  # SNP-thinning improves ascertainment of population structure with PCA
  # the following automatic procedure prunes and removes long-range LD regions (Privé, Aschard, and Blum 2017)
  # the snp_autoSVD function combines the following:
  #    + clumping of SNPs - consistent with the PLINK pruning method, whereby SNPs are ranked by MAF
  #    + but clumping keeps the most sig. SNP per region, rather than PLINK pruning which keeps the SNP with the highest MAF from a correlated pair
  #    + the authors reccomend use of clumping over pruning because pruning can result in empty regions of the genome
  #    + they note that in reality, the clumping and pruning algorithms produce similar sets of SNPs
  #    + this function also detects long-range LD which would cause SNPs to have consecutive loadings across PCs
  #    + SNPs in long-range LD are discarded before computing PCs
  # Arguments:
  #    + thr.r2 = Threshold over the squared correlation between two SNPs. Default is 0.2. Use NA if you want to skip the clumping step., in the paper, 0.05 is very stringent
  #    + k = number of PCs to calculate
  # The authors reccomend:
  #    + Always use both pruning (clumping) and removing of long-range LD regions when computing the PCs, as recommended by Abdellaoui et al. (2013).
  #    + To check that the results are capturing population structure, you should plot PCA scores.
  #    + To check that PCs are not capturing LD, you should check PCA loadings.
  svd0 <- snp_autoSVD(G, CHR, POS, thr.r2 = 0.2, k=10)
  
  print("Saving image")
  
  save.image(file = paste0(data.out, cohort, "_bigsnpr_data.RData"))
  
  ## 5. Export gPCs
  
  ## 5a. Get sample IDs
  sample_ids = obj.bigSNP$fam %>% 
    dplyr::pull(sample.ID)
  
  ## 5b. Add sample ids to sample-PC df and write out
  svd0$u %>%
    as.data.frame() %>% 
    `rownames<-`(sample_ids) %>% 
    readr::write_csv(x=., file=paste0(data.out, cohort, "_svd.obj$u_10_gPCs.csv"))

  print("----------------------")
  
}



