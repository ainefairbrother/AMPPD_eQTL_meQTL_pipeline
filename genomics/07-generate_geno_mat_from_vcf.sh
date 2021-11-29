#!/usr/bin/env bash

module load apps/vcftools

cd /scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files

echo $(ls *all_chrs.vcf.gz)

# this needs to be set to 3000 to prevent vcftools failing due to inability to open x number of temp files 
ulimit -n 3000

# set f_prefix to the file pattern suffix you want to recognise - i.e. the loop runs with all files with the suffix f_suffix
f_suffix="all_chrs-maf0.05-filtered.vcf.recode.vcf.gz"

for file in $(ls *${f_suffix});
do
fname=${file%.vcf.gz};
echo "Generating $fname ...";
# generate 012 matrix where: -1 represents missing data, 0 means 0 non-reference alleles, 1 means 1 non-ref allele, 2 means 2 non-ref alleles
vcftools --gzvcf $file --012 --out "${fname}_geno_matrix";
echo "$fname generated.";
done


