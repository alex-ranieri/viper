#!/usr/bin/env python
# coding: utf-8

# Importing libraries
import pandas as pd
import sys

# Reading input
statistics_file = sys.argv[1]
folder_name = sys.argv[2]

statistics_df = pd.read_csv(statistics_file, sep = '\t')

statistics_df_useful = statistics_df[['Genome', 'Segments_Assembled',
                                      'Type', 'H_genotype', 'N_genotype', 'Subtype', 'Clade',
                                     'segment_1_Coverage', 'segment_2_Coverage',
                                     'segment_3_Coverage', 'segment_4_Coverage',
                                     'segment_5_Coverage', 'segment_6_Coverage',
                                     'segment_7_Coverage', 'segment_8_Coverage']]
segments_series = {}
for index, row in statistics_df_useful.iterrows():
    seg_string = ''
    for i in ['segment_1_Coverage','segment_2_Coverage','segment_3_Coverage',
              'segment_4_Coverage','segment_5_Coverage','segment_6_Coverage','segment_7_Coverage',
              'segment_8_Coverage']:
        if (('not_detected' not in str(row[i])) and ('could_not_be_assembled' not in str(row[i]) )):
            seg_id = i.split('_')[1]
            seg_string = seg_string + '|' + seg_id
        segments_series[index] = seg_string
statistics_df_useful['Segments'] = pd.Series(segments_series)
passed_QC = {}
for index, row in statistics_df_useful.iterrows():
    seg4_detect = row['segment_4_Coverage']

    if ((seg4_detect == 'not_detected') or (seg4_detect == 'could_not_be_assembled')):
        passed_QC[index] = 'F'
    else:
        passed_QC[index] = 'A'
statistics_df_useful['Passed QC?'] = pd.Series(passed_QC)

statistics_df_useful['Notes']  = ''
statistics_df_useful['CEVIVAS_ID'] = ''

#criando a coluna repeat, estudo e analysis
statistics_df_useful['Repeat'] = ''
statistics_df_useful['Study'] = ''
statistics_df_useful['Analysis'] = ''

cevivas_output = statistics_df_useful[['CEVIVAS_ID','Genome', 'Type', 'Subtype', 'Clade', 'Segments_Assembled',
                        'Segments', 'Passed QC?', 'Notes', 'Repeat', 'Study', 'Analysis' ,'segment_1_Coverage','segment_2_Coverage',
                         'segment_3_Coverage','segment_4_Coverage','segment_5_Coverage',
                         'segment_6_Coverage','segment_7_Coverage','segment_8_Coverage']]
output_file_name = 'CeVIVAS_FLU_' + folder_name + '.xlsx'
cevivas_output.to_excel(output_file_name, index = False)

