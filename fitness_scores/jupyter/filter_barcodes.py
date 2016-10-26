import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import fitness_lib

R1_CUTOFF = 0.343451896717
R2_CUTOFF = 0.362603141536
CONTROL_CUTOFF = 0.324675892613


def filter_barcodes(barcode_scores_no_nans, cutoff_score):
	barcode_scores_filtered = {}
	for index_tup in barcode_scores_no_nans:
		if barcode_scores_no_nans[index_tup][1] < cutoff_score:
			barcode_scores_filtered[index_tup] = barcode_scores_no_nans[index_tup]

	return barcode_scores_filtered



barcode_scores_r1_no_nan = pickle.load(open('./input_output_files/output/r1_full_data_barcodes_no_nans.pkl', 'rb'))
barcode_scores_r2_no_nan = pickle.load(open('./input_output_files/output/r2_full_data_barcodes_no_nans.pkl', 'rb'))
barcode_scores_filtered_r1 = filter_barcodes(barcode_scores_r1_no_nan, R1_CUTOFF)
barcode_scores_filtered_r2 = filter_barcodes(barcode_scores_r2_no_nan, R2_CUTOFF)

# print barcode_scores_filtered_r1
pickle.dump(barcode_scores_filtered_r1, open('./input_output_files/output/r1_full_data_barcodes_filtered.pkl', 'wb'))
pickle.dump(barcode_scores_filtered_r2, open('./input_output_files/output/r2_full_data_barcodes_filtered.pkl', 'wb'))

