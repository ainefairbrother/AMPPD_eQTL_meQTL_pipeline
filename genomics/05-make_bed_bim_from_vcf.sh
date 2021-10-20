#!/bin/sh

#SBATCH --job-name=vcf-to-bed
#SBATCH --output=/scratch/groups/hodgkinsonlab/aine/data/per-sample-vcf/vcf-to-bed.log
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=300000

# date: 11/10/21
# purpose of script: to convert vcfs to bed files

module load apps/plink2

# go to relevant dir
cd /scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files

plink2 --vcf PP_all_chrs.vcf.gz --make-bed --out PP_all_chrs
plink2 --vcf PD_all_chrs.vcf.gz --make-bed --out PD_all_chrs
plink2 --vcf BF_all_chrs.vcf.gz --make-bed --out BF_all_chrs
