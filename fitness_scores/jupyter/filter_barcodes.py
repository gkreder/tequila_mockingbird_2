import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import fitness_lib

# -----------------------------------------------------------------------------------------------------------------------
# Filter the barcode scores
# -----------------------------------------------------------------------------------------------------------------------

# R1_CUTOFF_10 = 0.343451896717
# R2_CUTOFF_10 = 0.362603141536
# CONTROL_CUTOFF_10 = 0.324675892613

R1_CUTOFF_20 = 0.253472503491
R2_CUTOFF_20 = 0.28227116017
CONTROL_CUTOFF_20 = 0.28227116017


# def filter_barcodes(barcode_scores_no_nans, cutoff_score):
# 	barcode_scores_filtered = {}
# 	for index_tup in barcode_scores_no_nans:
# 		if barcode_scores_no_nans[index_tup][1] < cutoff_score:
# 			barcode_scores_filtered[index_tup] = barcode_scores_no_nans[index_tup]

# 	return barcode_scores_filtered



# barcode_scores_r1_no_nan = pickle.load(open('./input_output_files/output/r1_full_data_barcodes_no_nans.pkl', 'rb'))
# barcode_scores_r2_no_nan = pickle.load(open('./input_output_files/output/r2_full_data_barcodes_no_nans.pkl', 'rb'))
# barcode_scores_control_no_nan = pickle.load(open('./input_output_files/output/control_full_data_barcodes_no_nans.pkl', 'rb'))
# barcode_scores_filtered_r1 = filter_barcodes(barcode_scores_r1_no_nan, R1_CUTOFF_20)
# barcode_scores_filtered_r2 = filter_barcodes(barcode_scores_r2_no_nan, R2_CUTOFF_20)
# barcode_scores_filtered_control = filter_barcodes(barcode_scores_control_no_nan, R2_CUTOFF_20)

# # # print barcode_scores_filtered_r1
# pickle.dump(barcode_scores_filtered_r1, open('./input_output_files/output/r1_full_data_barcodes_filtered_20.pkl', 'wb'))
# pickle.dump(barcode_scores_filtered_r2, open('./input_output_files/output/r2_full_data_barcodes_filtered_20.pkl', 'wb'))
# pickle.dump(barcode_scores_filtered_control, open('./input_output_files/output/control_full_data_barcodes_filtered_20.pkl', 'wb'))

# -----------------------------------------------------------------------------------------------------------------------
# Get filtered AA scores and write to file
# -----------------------------------------------------------------------------------------------------------------------

barcode_scores_filtered_r1 = pickle.load(open('./input_output_files/output/r1_full_data_barcodes_filtered_20.pkl', 'rb'))
barcode_scores_filtered_r2 = pickle.load(open('./input_output_files/output/r2_full_data_barcodes_filtered_20.pkl', 'rb'))
barcode_scores_filtered_control = pickle.load(open('./input_output_files/output/control_full_data_barcodes_filtered_20.pkl', 'rb'))


AA_scores_filtered_r1, thrown_out_dictionary_reads_r1 = fitness_lib.map_barcode_to_AA(barcode_scores_filtered_r1)
AA_scores_filtered_r2, thrown_out_dictionary_reads_r2 = fitness_lib.map_barcode_to_AA(barcode_scores_filtered_r2)
AA_scores_filtered_control, thrown_out_dictionary_reads_control = fitness_lib.map_barcode_to_AA(barcode_scores_filtered_control)

pickle.dump(AA_scores_filtered_r1, open('./input_output_files/output/r1_filtered_AA_scores_20' , 'wb'))
pickle.dump(AA_scores_filtered_r2, open('./input_output_files/output/r2_filtered_AA_scores_20' , 'wb'))
pickle.dump(AA_scores_filtered_control, open('./input_output_files/output/control_filtered_AA_scores_20' , 'wb'))

# -----------------------------------------------------------------------------------------------------------------------
# Try out heatmap
# -----------------------------------------------------------------------------------------------------------------------


