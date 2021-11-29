#!/bin/sh

#SBATCH --job-name=filter-vcfs
#SBATCH --output=/scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files/vcf-filter-log.log
#SBATCH --time=0-48:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=150000

# date created: 25/11/21
# last updated: 25/11/21
# purpose of script: to filter a directory of vcf files and save smaller filtered files in the same directory
# filter split cohort, rejoined chr files for MAF<5%

# load modules
#module load apps/vcftools
#module load apps/samtools

# set working dir 
wd="/scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files"

# go to working dir 
cd ${wd}

# set maf
maf=0.05

echo "Running loop on the following files in" ${wd}
echo $(ls *all_chrs.vcf.gz)

for i in $(ls *all_chrs.vcf.gz);
do

filename="${i%%.*}"
outname=${filename}-maf${maf}-filtered.vcf
echo "Current file: " $i "-->" $outname

echo "1. Filtering file"

vcftools --gzvcf ${i} --maf $maf --recode --recode-INFO-all --out $outname

echo "2. Compressing"
bgzip -c $outname > $outname.gz;

echo "3. Indexing"
tabix -p vcf $outname.gz;

done


