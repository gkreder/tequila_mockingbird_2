import fitness_lib
import pickle
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from operator import itemgetter

AMINOTONUMBER_DATA = pickle.load(open('./input_output_files/input/aminotonumber.pkl', 'rb'))
scriptname, replicate = sys.argv
file_name = './input_output_files/output/' + replicate + '_full_data.pkl'


def plot_hmap(data, row_labels, column_labels):
    # sns.set(font_scale=0.7)
    sns.set(font_scale=1.0)

    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(data, ax=ax,
                     cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"},
                     xticklabels=row_labels,
                     yticklabels=column_labels,
                     linewidths = 0.5)
    
    # mpl.rcParams.update({'font.size': 22})
    plt.sca(ax)
    plt.yticks(rotation=0)
    plt.xticks(rotation=0)
    plt.title(replicate + ' Fitness Scores')
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





# (AA_scores, total_reads, thrown_out_N_reads, thrown_out_dictionary_reads) = pickle.load(open(file_name, 'rb'))
# plot_scores(AA_scores)