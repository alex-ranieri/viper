#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import pandas as pd
import sys
import numpy as np

# Getting files from sys.argv
stats_file = sys.argv[1]
nextclade_file = sys.argv[2]

# Reading files 
stats = pd.read_csv(stats_file, sep='\t')
nextclade = pd.read_csv(nextclade_file, sep='\t')

# Processing
nextclade = nextclade.rename(columns = {'seqName':'Genome'}) # Renaming column in nextclade output
#nextclade['Genome'] = ["Genoma_" + str(x) + ".fasta" for x in nextclade["Genome"]] # Adding "Genoma_" prefix and ".fasta" sufix 
nextclade_useful = nextclade[['Genome','clade','qc.overallStatus', 'deletions','insertions','aaSubstitutions', 'aaDeletions', 'totalFrameShifts', 'frameShifts', 'substitutions']] # Selecting only necessary information
#df = pd.merge(stats, nextclade_useful, on=['Genome']) # Merging data
df = pd.merge(stats, nextclade_useful,on=['Genome'], how="outer")
#df = df[df['clade'].notna()] # Dropping data without clade information
#df['Passed QC?'] = np.where(df['Npos_Depth>=10'] < 20000, 'N', 'S')
df['Passed QC?'] = np.where(df['Npos_Depth>=10'] < 20000, 'N', np.where(df['Coverage'] < 85, 'N', np.where(df['totalFrameShifts'] > 0, 'M', np.where(df['SNPs_count'] >= 5, 'R', 'S'))))

output_fileName = stats_file[0:-11] + ".xlsx"
df.to_excel(output_fileName, index=False) # Writing results to file 
