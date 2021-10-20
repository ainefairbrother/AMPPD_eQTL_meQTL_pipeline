#!/bin/sh
#$ -cwd
#$ -V

cd /scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files

for i in $(ls *.vcf);
do

echo ${i};

echo "Compressing" ${i} "-->" ${i}.gz; 
bgzip -c $i > ${i}.gz;

echo "Indexing" ${filename}.gz; 
tabix -p vcf ${i}.gz;

done

