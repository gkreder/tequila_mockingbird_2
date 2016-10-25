import os, sys
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import pickle
import numpy as np
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import warnings
import fitness_lib
import plot_scores


scriptname, replicate = sys.argv
def unpack_AA_scores(filename):
	(AA_scores, total_reads, thrown_out_N_reads, thrown_out_dictionary_reads) = pickle.load(open(filename, 'rb'))
	return AA_scores

def unpack_barcode_scores(filename):
	barcode_scores = pickle.load(open(filename, 'rb'))
	return barcode_scores
	# print barcode_scores

def filter_nans(barcode_scores):
	filtered_barcode_scores = {}
	for barcode in barcode_scores:
		if not np.isnan(barcode_scores[barcode][0]):
			filtered_barcode_scores[barcode] = barcode_scores[barcode]

	print filtered_barcode_scores
	return filtered_barcode_scores

def plot_hist(filtered_barcode_scores, metric = 'score'):
	values = []
	if metric == 'score':
		index = 0
	else:
		index = 1

	for tup in list(filtered_barcode_scores.values()):
		values.append(tup[index])

	plt.hist(values)
	plt.show()





replicate = 'r1'
fname_AA = './input_output_files/output/' + replicate + '_full_data.pkl'
fname_barcode = './input_output_files/output/' + replicate + '_full_data_barcodes.pkl'
# unpack_AA_scores(fname)
r1_barcode_scores = unpack_barcode_scores(fname_barcode)
r1_filtered_barcode_scores = filter_nans(r1_barcode_scores)
# AA_scores_filtered, thrown_out_dictionary_reads_filtered = fitness_lib.map_barcode_to_AA(r1_filtered_barcode_scores)
plot_hist(r1_filtered_barcode_scores)





