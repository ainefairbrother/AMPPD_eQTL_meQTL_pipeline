#!/bin/sh

#SBATCH --job-name=split-vcf
#SBATCH --output=/scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files/split-vcfs-by-sample.log
#SBATCH --time=2-00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=50000

# date: 14/10/21
# purpose of script: to split per-chr vcf files into cohorts
# run the script generate_cohort_sampleid_lists.sh first to generate the sample lists required to split the files by cohort 

conda activate sambcfenv

# go to relevant dir
cd /scratch/groups/hodgkinsonlab/aine/data/amppd_vcf_files

echo $(pwd)

declare -a cohorts=("PP" "PD" "BF")

echo $cohorts

echo "Starting analysis loop ..."

# list of files to perform sample filtration on 
files=$(ls *maf0.01-missing0.01-hwe0.001-filtered*)

"for file in ${files[@]};
do
      echo "file: "$file;

      for cohort in ${cohorts[@]};
      do
            echo "cohort: "$cohort";
            sample_list="./sample_lists_cohort_chr/${file%.*.*}-${cohort}.txt";   
	    outfile_name="${file%.*.*}-${cohort}.vcf";
            echo "sample list: "$sample_list;
            echo "output file name: "$outfile_name;
	    bcftools view -S $sample_list -Oz -o $outfile_name $file; 
	done
done
