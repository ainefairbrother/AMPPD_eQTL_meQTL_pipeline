#!/bin/sh

#SBATCH --job-name=run-peer
#SBATCH --output=/scratch/groups/hodgkinsonlab/aine/data/amppd_log10mednorm/run-peer.log
#SBATCH --time=0-24:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20480

# date: 28/09/21
# purpose of script: to run peer

fname="PP_ctrl"
indir="/scratch/groups/hodgkinsonlab/aine/data/amppd_log10mednorm/"

cd $indir
mkdir $fname
cd $fname

/scratch/groups/hodgkinsonlab/eMedLab/away/ahodgkinson/Programs/peer/software/bin/peertool -f ${indir}${fname}_log10_mediannorm_TPM_wrangled.csv -n 50 -i 1000 --has_header --has_rownames -o .


