#!/bin/sh
#$ -cwd
#$ -V

# this is a script to check filtered vcf files for validity 
# it will produce a report for each file detailing:
# 1. number of samples
# 2. number of lines in this file vs. in original (pre-filtered) file 
# 3. final line of the file 
# 4. column names 

# need to run this in: 
# conda activate sambcfenv

cd /scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files

for i in $(ls *recode*);
do

	touch ${filename}.validation.report.txt

	echo ${i}
	filename="${i%%-*}"
	echo ${filename}

	# 1. 

	echo "File: ${filename}" >> ${filename}.validation.report.txt

	echo "1. Number of samples" >> ${filename}.validation.report.txt

	bcftools query -l $i | wc -l >> ${filename}.validation.report.txt

	# 2. 

	#echo "2. Count number of lines" >> ${filename}.validation.report.txt

	#echo "i. in the original file" >> ${filename}.validation.report.txt
	#bcftools query -l ${filename}.vcf.gz | sort | uniq | wc -l >> ${filename}.validation.report.txt
	#bcftools view ${filename}.vcf.gz | wc -l | cut -d' ' -f1  >> ${filename}.validation.report.txt

	#echo "ii. in the filtered file" >> ${filename}.validation.report.txt
	#bcftools query -l $i | sort | uniq | wc -l >> ${filename}.validation.report.txt
	#bcftools view ${i} | wc -l | cut -d' ' -f1 >> ${filename}.validation.report.txt

	# 3. 

	# 4.

	echo "4. get %CHROM  %POS  %ID  %REF  %ALT{0}  %QUAL cols and print 4 lines of each" >> ${filename}.validation.report.txt
	bcftools query -f '%CHROM  %POS  %ID  %REF  %ALT{0}  %QUAL\n' ${i} | head -n 4 >> ${filename}.validation.report.txt
done

