#!/usr/bin/env python3

# To run this script, methylation matrices (1 file per sample, 1 row per file) were transferred from the Terra platform using gsutils to the relevant server
# The matrices were generated from AMP-PD BAM files using the following workflow: https://portal.firecloud.org/?return=terra#methods/aine_fb_ucl/calculate_mismatch_proportion_workflow/14

import csv
import pandas as pd
import os
import re
import glob

path = r'/scratch/groups/hodgkinsonlab/aine/data/amppd_meth_matrices/'

# get list of files in dir
all_files = [path+file for file in os.listdir(path) if 'MTmethmatrix' in file]

print('1. Making generator object')
# make generator object, reading in each file
df_from_each_file = (pd.read_csv(f, sep='\t', encoding="utf-8", engine='python', index_col=0, header=0) for f in all_files)

print('2. Concatenating dfs')
# concatenate
concatenated_df = pd.concat(df_from_each_file, ignore_index=True)

print('3. Adding sample_id col')
# insert sample id col and drop IID
concatenated_df.insert(loc=0, column='sample_id', value=[re.search('.+\/(.+)\.MTmethmatrix\.csv', f).group(1) for f in all_files])
concatenated_df.drop(['IID'], axis=1, inplace=True)

print('4. Saving...')
# save concatenated df 
concatenated_df.to_csv(path+'amppd_methylation_matrix.csv', index=True, header=True)
