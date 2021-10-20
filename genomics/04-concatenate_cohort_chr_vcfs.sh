#!/bin/sh
# date: 11/10/21

# purpose of script: to concatenate chr-cohort vcf files into all_chr cohort vcf files

conda activate sambcfenv

declare -a cohorts=("PP" "PD" "BF")

for cohort in ${cohorts[@]}; 
do 
	# identify files with chrN and cohortX, then combine for each cohort into one vcf
	cohort_files=$(ls -1v chr*${cohort}.vcf); 
	echo $cohort; 
	echo $cohort_files; 
	bcftools concat --threads 4 -o ${cohort}_all_chrs.vcf.gz $cohort_files; 
done


