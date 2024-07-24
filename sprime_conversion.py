import pandas as pd
import numpy as np
import os
import glob
from pathlib import PureWindowsPath

##read in all sprime files and create one big df - this directory has all sprime files
path = r'C:\\Users\\SciFunk\\Downloads\\Working\\Biostats\\sprimeResults\\'              
all_files = glob.glob(os.path.join(path, "*"))

df_from_each_file = (pd.read_csv(f, sep='\t').assign(pop=PureWindowsPath(f).name.split('.')[0]) for f in all_files)
chr_df = pd.concat(df_from_each_file, ignore_index=True)

##set a new column for the archaic allele so we can drop 'REF', 'ALT', 'ALLELE'
chr_df['arc_allele'] = np.where((chr_df['ALLELE'] == 1), chr_df['ALT'], chr_df['REF'])
chr_df.drop(['REF', 'ALT', 'ALLELE', 'SEGMENT', 'SCORE'], axis = 1, inplace = True)

allArchaic = chr_df.loc[(chr_df['NMATCH'] == "match") | (chr_df['DMATCH'] == "match")] 

##create a properly formatted bedfile for samtools
allArchaic.drop_duplicates(['ID'], inplace = True)
allArchaic['chrom'] = allArchaic['CHROM']
allArchaic['chromStart'] = allArchaic['POS']
allArchaic['chromEnd'] = allArchaic['POS'] + 1

allArchaic[['chrom', 'chromStart', 'chromEnd']]
allArchaic[['chrom', 'chromStart', 'chromEnd']].to_csv('sPrime_sites.bed', sep="\t")
