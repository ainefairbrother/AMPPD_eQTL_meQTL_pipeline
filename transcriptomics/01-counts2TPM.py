
# allows fn to be run from the cmd line
import fire

def counts2TPM(file_dir="", gene_length_dir= "", pattern=".csv", filesep=",", outdir=None, outlabel=""):
  
   # """
  # file_dir: directory containing the rpkms
  # pattern: pattern to search for, default is all .csv files
  # outdir: where should the files go, default is file_dir
  # filesep: file separating char, default is comma
  # outlabel: add label to output files
  # """
  
  import os 
  import pandas as pd
  
  if outdir == None:
    outdir = file_dir
  
  # get all file paths
  file_paths = [file_dir+file for file in os.listdir(file_dir) if pattern in file]
  file_names = [file for file in os.listdir(file_dir) if pattern in file]
  print(file_names)
  
  for i in range(0,len(file_paths)):
    
    print(file_names[i])
    
    # import df
    df = pd.read_csv(open(file_paths[i]), sep=filesep, encoding="utf-8", engine='python', index_col=0, header=0) # genes cols, samples rows
    
    # genes to rows 
    if len(df.index.values) < len(df.columns.values):
      df = df.T
    print("counts shape = {}".format(df.shape))
    
    # import gene lengths
    gene_lengths = pd.read_csv(open(gene_length_dir), encoding="utf-8", engine='python', index_col=0, header=0)
    gene_lengths.reset_index(inplace=True)
    gene_lengths.columns = ['gene_id', 'width']
    
    print("-----------------------------")
    print("df")
    print(df.head())
    print(df.index.values[1:5])
    print("gene_lengths")
    print(gene_lengths.head())
    print("-----------------------------")
    
    # filter gene lengths df for genes in counts file
    gene_lengths_filtered = gene_lengths[gene_lengths['gene_id'].isin(df.index.values.tolist())]
    print("filtered gene length shape = {}".format(gene_lengths.shape))
    gene_length_dict = dict(zip(gene_lengths_filtered['gene_id'], gene_lengths_filtered['width']))
    
    # join tables
    df["width"] = df.index.to_series().map(gene_length_dict)
    
    # 1. RPK (reads per kilobase) = divide the read counts by the length of each gene in kb and include only rows w/ relevant sample name label
    label_starting_chars = df.columns.tolist()[0][0:2]
    RPK = df.filter(regex=label_starting_chars).div(df['width'], axis=0)
    
    # drop width col
    if 'width' in df.columns.values:
      del df['width']
    
    # 2. PMSF = count all the RPK values in a sample and divide this number by 1m, this is the "per million scaling factor"
    per_million_scaling_factor = RPK.sum(axis=0)/1000000
    
    # 3. Calculating TPM from RPK
    TPM = RPK / per_million_scaling_factor
    
    # genes to cols
    #TPM = TPM.T
    
    print(TPM.head())

    # export filtered file
    os.chdir(outdir)
    TPM.to_csv(outlabel + "_TPM.csv", index=True, header=True)

# allows fn to be run from the cmd line
if __name__ == '__main__':
  fire.Fire(counts2TPM)
    
    
    
    
    
    
    
    
    
