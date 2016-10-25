import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import pickle
import numpy as np
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
# import fastq scrape

time_points = [
    't0_r1',
    't1_r1',
    't2_r1',
    't0_r2',
    't1_r2',
    't2_r2',
    't0_control',
    't1_control',
    't2_control']

ID_timepoint_dict = {
    'ATTCCG' : 't0_control',
    'AGCTAG' : 't1_control',
    'GTATAG' : 't2_control',
    'ATCAGT' : 't0_r1',
    'GCTCAT' : 't1_r1',
    'AGGAAT' : 't2_r1',
    'CTTTTG' : 't0_r2',
    'TAGTTG' : 't1_r2',
    'CCGGTG' : 't2_r2'
}

timepoint_generation_dict = {'t0_r1' : 0.0,
                             't1_r1' : 1.91,
                             't2_r1': 3.328354364,
                             't0_r2' : 0.0,
                             't1_r2' : 1.35,
                             't2_r2' : 2.96,
                             't0_control' : 0.0,
                             't1_control' : 2.022367813,
                             't2_control' : 3.50041511
                            }


allele_data = pickle.load(open('allele_dic_with_WT.pkl', 'rb'))
aminotonumber_data = pickle.load(open('aminotonumber.pkl', 'rb'))
translate_data_RNA = pickle.load(open('translate.pkl', 'rb'))
translate_data_DNA = {}
for codon in translate_data_RNA:
    translate_data_DNA[codon.replace('U', 'T')] = translate_data_RNA[codon]
translate_data_DNA['WT'] = 'WT'

def compare_string_lists(l1, l2):
#     print(l1, l2)
    for index, element in enumerate(l1):
        if len(l1) != len(l2):
            return False
        if element != l2[index]:
            return False
    return True



# 
# Initial Data Scrape from Fastq files
# 
df = pd.DataFrame()
ID = []
read = []
quality_score = []
barcodes = []
codon = []
loc = []
time_point_tag = []

def get_filename(time_point):
    # filename_short = '../../' + time_point + '_short.fastq'
    # return filename_short
    filename_long =  time_point + '_index.fastq'
    return filename_long
    
for time_point in time_points:
    filename = get_filename(time_point)
    file = open(filename)
    for line_number, line in enumerate(file):
        if line_number % 4 == 0:
            ID.append(line.strip().split(':')[-1])
        elif line_number % 4 == 1:
            read.append(line.strip())
            barcode = line.strip()[0:18]
            barcodes.append(barcode)
            barcode_lookup = str(Seq(barcode).reverse_complement())
            if barcode_lookup in allele_data:
                codon.append(allele_data[barcode_lookup][1])
                loc.append(allele_data[barcode_lookup][0])
            else:
                codon.append('NONE')
                loc.append('NONE')

        elif line_number % 4 == 3:
            quality_score.append(line.strip())
            time_point_tag.append(time_point)
df['ID'] = ID
df['read'] = read
df['time_point'] = time_point_tag
df['codon'] = codon
df['loc'] = loc



# 
# Get normalized barcode counts per timestep and label with generation, timepoint, and replicate
# 

df_in_dict = df[~(df['codon'] == 'NONE')].copy()
# df_in_dict.groupby('time_point').codon.nunique()
grouped_counts = df_in_dict.groupby('time_point').read.value_counts(normalize = True)
keys = grouped_counts.keys()
def label_row_counts(row):
    temp_key = (row['time_point'], row['read'])
    for key in keys:
        if temp_key == key:
            return grouped_counts[key]
    return 'NO COUNT FOUND'

def label_row_generations(row):
    return timepoint_generation_dict[row['time_point']]

def label_row_replicate(row):
    return row['time_point'].split('_')[1]
def label_row_time(row):
    return row['time_point'].split('_')[0]

df_in_dict['count'] = df_in_dict.apply(lambda row: label_row_counts(row), axis=1)
AA_list = [translate_data_DNA[codon] for codon in df_in_dict['codon'].values.tolist()]
df_in_dict['AA'] = AA_list
# df_in_dict[df_in_dict['time_point'] == 't2_r2']
df_total = df_in_dict.copy()
df_total['generation'] = df_total.apply(lambda row: label_row_generations(row), axis=1)
df_total['time'] = df_total.apply(lambda row: label_row_time(row), axis=1)
df_total['replicate'] = df_total.apply(lambda row: label_row_replicate(row), axis=1)



# 
# Calculate Slopes (can change our filtering of slope values here)
# 
df_slopes = df_total.drop_duplicates().copy()
barcode_groups = df_slopes.groupby(['read', 'replicate'])
groups = barcode_groups.groups
for group in barcode_groups.groups:
    replicate = group[1]
    temp_group = barcode_groups.get_group(group)
    temp_x = temp_group['generation'].values
    temp_y = temp_group['count'].values
    timepoints_present = temp_group['time'].values
    if compare_string_lists(timepoints_present, ['t1','t2']):
        temp_x = [0.0] + temp_x
        temp_y = [0.0] + temp_y
    elif compare_string_lists(timepoints_present, ['t0']):
        generation_1 = timepoint_generation_dict['t1_' + replicate]
        generation_2 = timepoint_generation_dict['t2_' + replicate]
        temp_x = temp_x + [generation_1, generation_2]
        temp_y = temp_y + [0.0, 0.0]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(temp_x,temp_y)
    indices_temp = temp_group.index.values
    for index in indices_temp:
        df_slopes.set_value(index, 'slope', slope)
        df_slopes.set_value(index, 'std_err', std_err)



# 
# Calculate the average slope values
# 
df_barcode_scores = df_slopes.drop_duplicates(subset = ['read', 'replicate']).copy()
df_AA_avg = pd.DataFrame(df_barcode_scores.groupby(['loc', 'AA', 'replicate']).mean().copy())
pickle.dump(df_AA_avg, open('df_AA_avg.pkl', 'wb'))

# keys_temp = df_AA_avg['slope'].keys()
# for key in keys_temp:
#     print(key)


# 
# Separate out into replicate DFs
# 
def subtract_wildtype(df_replicate):
    WT_slope = df_replicate[df_replicate['AA'] == 'WT']['slope']
    normalized_slopes = df_replicate.apply(lambda row : row['slope'] - WT_slope, axis = 1)
    df_replicate['slope_normalized'] = normalized_slopes
    df_return = df_replicate.copy()
    return df_return

df_r1 = df_AA_avg.reset_index().copy()
df_r1 = pd.DataFrame(df_r1[df_r1['replicate'] == 'r1'])
df_r1 = subtract_wildtype(df_r1)
pickle.dump(df_r1, open('df_r1.pkl', 'wb'))

df_r2 = df_AA_avg.reset_index().copy()
df_r2 = pd.DataFrame(df_r2[df_r2['replicate'] == 'r2'])
df_r2 = subtract_wildtype(df_r2)
pickle.dump(df_r2, open('df_r2.pkl', 'wb'))

df_control = df_AA_avg.reset_index().copy()
df_control = pd.DataFrame(df_control[df_control['replicate'] == 'control'])
df_control = subtract_wildtype(df_control)
pickle.dump(df_control, open('df_control.pkl', 'wb'))





def plot_df_replicate(df, filename, normalized = True):
    df_working = df.copy()
    df_working = df_working.drop('std_err', 1)
    df_working = df_working.drop('count', 1)
    df_working = df_working.drop('generation', 1)
    df_working = df_working.drop('replicate', 1)
    if normalized:
        df_working = df_working.drop('slope', 1)
        df_working = df_working.pivot(index='AA', columns='loc', values='slope_normalized')
    else:
        df_working = df_working.drop('slope_normalized', 1)
        df_working = df_working.pivot(index='AA', columns='loc', values='slope')
    pickle.dump(df_working, open(filename, 'wb'))
#     sns.heatmap(df_working)
#     plt.show()

plot_df_replicate(df_control, 'control_heat_df.pkl')
plot_df_replicate(df_r1, 'r1_heat_df.pkl')
plot_df_replicate(df_r2, 'r2_heat_df.pkl')