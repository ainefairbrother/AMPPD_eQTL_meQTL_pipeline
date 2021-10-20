#!/bin/sh

cd /scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files

declare -a cohorts=("PP" "PD" "BF") 
	
for file in $(ls chr*-filtered.vcf.gz); 
do
	for cohort in ${cohorts[@]}; 
		do bcftools query -l $file | grep "${cohort}-" > ./sample_lists_cohort_chr/${file%.*.*}-${cohort}.txt;
	done
done



