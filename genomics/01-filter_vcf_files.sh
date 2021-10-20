#!/bin/sh

#SBATCH --job-name=filter-vcfs
#SBATCH --output=/scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files/vcf-filter-log.log
#SBATCH --time=0-168:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=300000

# date: 05/08/21
# purpose of script: to filter a directory of vcf files and save smaller filtered files in the same directory
# vm needed to run/ prereqs: use `conda activate st_new`
# CMD line args:
# $1 corresponds to the dir in which the vcf files are stored

module load apps/vcftools

wd="/scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files"

cd ${wd}

touch TESTFILE.txt

echo "Running loop on the following files in" ${wd}
echo $(ls -lah *.vcf.gz)

for i in $(ls *.vcf.gz);
do

filename="${i%%.*}"
echo "Current file: " ${i} "-->" ${filename}-maf0.01-missing0.01-hwe0.001-filtered.vcf.gz

echo "1. Filtering file"
#vcftools --gzvcf ${i} --maf 0.01 --max-missing 0.01 -–hwe 0.001 --recode --recode-INFO-all -std-out | bgzip -c ${filename}-maf0.01-missing0.01-hwe0.001-filtered.vcf > ${filename}-maf0.01-missing0.01-hwe0.001-filtered.vcf.gz

vcftools --gzvcf ${i} --maf 0.01 --max-missing 0.01 --hwe 0.001 --recode --recode-INFO-all --out ${filename}-maf0.01-missing0.01-hwe0.001-filtered.vcf

done
