#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import pandas as pd
import sys

# Reading input
statistics_file = sys.argv[1]
folder_name = sys.argv[2]

#carregar os arquivos
statistics_df = pd.read_csv(statistics_file, sep= '\t')

#estruturar 
statistics_df_useful = statistics_df[['Genome', 'Serotype', 'Genotype', 'SNPs', 'Coverage', 'E_Coverage']]

#gerar o passed qc
passed_QC = {}
for index, row in statistics_df_useful.iterrows():
    cov = row['Coverage']
    ecov = row['E_Coverage']
    if (cov) >= 85.0:
        passed_QC[index] = 'A'
    elif (ecov) >= 55.0:
        passed_QC[index] = 'E'
    else:
        passed_QC[index] = 'R'
statistics_df_useful['Passed QC?'] = pd.Series(passed_QC)

#criando a coluna notes 
statistics_df_useful['Notes'] = ''

#criando a coluna cevivas ID
statistics_df_useful['CEVIVAS_ID'] = ''

#criando a coluna repeat, estudo e analysis
statistics_df_useful['Repeat'] = ''
statistics_df_useful['Estudo'] = ''
statistics_df_useful['Analysis'] = ''

#selecionando as coluna e colando-as em ordem 
statistics_df_useful = statistics_df_useful[['CEVIVAS_ID','Genome','Serotype','Genotype','Coverage', 'E_Coverage', 'Passed QC?','SNPs',
                                 'Notes','Repeat','Estudo','Analysis']]

#criando o arquivo final 
output_file_name = 'CeVIVAS_DENV_' + folder_name + '.csv'
statistics_df_useful.to_csv(output_file_name, index = False)
