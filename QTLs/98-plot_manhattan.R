#!/usr/bin/env Rscript

# generating manhattan plots from wrangled matrixeqtl output files

#### -------------- prepare env -------------- ####

# conda activate r4-base
# cd /home/abrowne/projects/amppd_analysis/pipeline/eqtls_mqtls
# Rscript 98-plot_manhattan.R

# to set up an env to run this script, do:
# conda install -c conda-forge r-stringr
# conda install -c conda-forge r-tidyr
# conda install -c conda-forge r-tibble
# conda install -c conda-forge r-ggplot2
# conda install -c conda-forge r-dplyr

#### -------------- load libs -------------- ####
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)

#### -------------- define functions -------------- ####

plot.manhattan.from.file = function(matrixeqtl.wrangled.file, print.or.save="save"){
  
  parse.filename = stringr::str_match(string=matrixeqtl.wrangled.file, pattern="\\/.+\\/(P.+)_(.+)_pheno=(.+)_maf=0\\.05_")
  cohort=parse.filename[,2]
  diagnosis=parse.filename[,3]
  phenotype=parse.filename[,4]
  
  if(file.exists(matrixeqtl.wrangled.file)==FALSE){
    print(paste("file", matrixeqtl.wrangled.file, "does not exist"))}
  else if(file.exists(matrixeqtl.wrangled.file)==TRUE){
    
    # filter the dataset for the input args
    gwasResults = vroom::vroom(matrixeqtl.wrangled.file) %>% 
      # replace X and Y chr labels with numbers
      dplyr::mutate(snp.chr = replace(snp.chr, snp.chr == "X", 23)) %>% 
      dplyr::mutate(snp.chr = replace(snp.chr, snp.chr == "Y", 24)) %>% 
      dplyr::mutate(snp.chr=as.numeric(snp.chr)) %>% 
      
      # rename for compatibility 
      dplyr::rename(
        # SNP=SNP,
        CHR=snp.chr,
        BP=snp.start,
        P=p.value # change to p.value for raw or FDR for corrected
      ) 
    
    if(dim(gwasResults)[1]!=0){
      
      don = gwasResults %>% 
        
        # Compute chromosome size
        dplyr::group_by(CHR) %>% 
        dplyr::summarise(chr_len=max(BP)) %>% 
        dplyr::mutate(chr_len=as.numeric(chr_len)) %>% 
        
        # Calculate cumulative position of each chromosome
        dplyr::mutate(tot=cumsum(chr_len)-chr_len) %>%
        dplyr::select(-chr_len) %>%
        
        # Add this info to the initial dataset
        dplyr::left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
        
        # Add a cumulative position of each SNP
        tidyr::drop_na(CHR) %>% 
        tidyr::drop_na(BP) %>% 
        dplyr::arrange(CHR, BP) %>% 
        dplyr::mutate( BPcum=BP+tot) %>% 
        
        # Add highlight and annotation information
        #mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
        mutate(is_annotate=ifelse(P<5e-8, "yes", "no")) 
      
      axisdf = don %>% 
        dplyr::group_by(CHR) %>% 
        dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
      
      plot = ggplot(don, aes(x=BPcum, y=-log10(P))) +
        
        # Show all points
        geom_point(aes(color=as.factor(CHR)), alpha=0.75, size=0.8) +
        scale_color_manual(values = rep(c("grey", "black"), 22)) +
        
        # custom X axis:
        scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
        scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +     # remove space between plot area and x axis
        
        # custom theme:
        theme_bw() +
        theme( 
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        ) + 
        geom_hline(yintercept = -log10(5e-8), colour="red", alpha=0.5)  + 
        xlab("Genomic position") +
        
        # custom theme:
        theme_bw(base_size = 8) +
        theme( 
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )
        ggtitle(paste(cohort, diagnosis, phenotype))
      
      # Add highlighted points
      #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
      
      # Add label using ggrepel to avoid overlapping
      #ggrepel::geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2)
      
      if(print.or.save=="print"){return(plot)}
      if(print.or.save=="save"){
        
        savePlot <- function(the.plot, out.name) {
          png(paste(out.name), width = 5, height = 3, units="in", res=200)
          print(the.plot)
          dev.off() 
        }
        
        savePlot(the.plot=plot, out.name=gsub("_wrangled.csv", "_manhattan.png", matrixeqtl.wrangled.file))
      }
    }
  }
}

lapply.helper = function(fn, path, pattern){
  
  # Run a function, fn, across files matching pattern, pattern, in directory, path.
  # This will produce whatever files fn is designed to output
  file.list = list.files(path=path, pattern=pattern, full.names=T)
  lapply(X=file.list, FUN=fn)
  
}

#### -------------- implement functions -------------- ####

lapply.helper(fn=plot.manhattan.from.file, 
              path="~/projects/amppd_analysis/data/MatrixEQTL_output/PP_PD_mega_analysis/", 
              pattern="_wrangled.csv")
