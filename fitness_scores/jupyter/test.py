import fitness_lib
import pickle
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from operator import itemgetter

AMINOTONUMBER_DATA = pickle.load(open('./input_output_files/input/aminotonumber.pkl', 'rb'))
replicate = 'control'
# file_name = './input_output_files/output/' + replicate + '_full_data.pkl'
file_name = './input_output_files/output/control_full_data.pkl'

# def plot_scores(AA_scores):
#     print(AA_scores)

def plot_hmap(data, row_labels, column_labels):
    sns.set(font_scale=0.7)

    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(data, ax=ax,
                     cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"},
                     xticklabels=row_labels,
                     yticklabels=column_labels,
                     linewidths = 0.5)
    plt.sca(ax)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90)
    plt.show()


def plot_scores(AA_scores):
    locs = []
    # AAs_keys = AMINOTONUMBER_DATA.keys().
    AA_labels = [tup[0] for tup in sorted(AMINOTONUMBER_DATA.items(), key=itemgetter(1))]
    scores = []

    # numpy.zeros((21,77))
    scores = np.empty((21,75))
    scores[:] = np.NAN
    
    for (loc, AA) in AA_scores:
        print loc, AA
        if AA != 'WT' and loc > 1 and loc < 77:
            scores[int(AMINOTONUMBER_DATA[AA]), int(loc - 2)] = AA_scores[(loc, AA)]['score']
            locs.append(loc)

    loc_labels = sorted(set(locs))
    plot_hmap(scores, loc_labels, AA_labels)



    # for (loc, AA) in AA_scores:
    #     locs.append(loc)
    #     AAs.append(AA)
    #     scores.append(AA_scores[(loc, AA)]['score'])
    # plotting_df = pd.DataFrame()
    
    # plotting_df['Location'] = locs
    # plotting_df['AA'] = AAs
    # plotting_df['score'] = scores
    # plotting_df = plotting_df[plotting_df['Location'] != 0]
    # plotting_df = plotting_df[plotting_df['AA'] != 'WT']
    # plotting_df = plotting_df.pivot(index='AA', columns='Location', values='score')
    # print plotting_df.sort()
    # # sns.heatmap(plotting_df)
    # # plt.show()






(AA_scores, total_reads, thrown_out_N_reads, thrown_out_dictionary_reads) = pickle.load(open(file_name, 'rb'))
plot_scores(AA_scores)
# plot_scores(AA_scores)