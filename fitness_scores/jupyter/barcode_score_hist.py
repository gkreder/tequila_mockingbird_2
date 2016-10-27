import pickle
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# # barcode_scores_r1 = pickle.load(open('./input_output_files/output/r1_full_data_barcodes.pkl', 'rb'))
# # barcode_scores_r2 = pickle.load(open('./input_output_files/output/r2_full_data_barcodes.pkl', 'rb'))
# barcode_scores_control = pickle.load(open('./input_output_files/output/control_full_data_barcodes.pkl', 'rb'))


# def filter_nans(barcode_scores):
# 	filtered_barcode_scores = {}
# 	for barcode in barcode_scores:
# 		if not np.isnan(barcode_scores[barcode][0]):
# 			filtered_barcode_scores[barcode] = barcode_scores[barcode]
# 	# print filtered_barcode_scores
# 	return filtered_barcode_scores

# # barcode_scores_r1_no_nan = filter_nans(barcode_scores_r1)
# # barcode_scores_r2_no_nan = filter_nans(barcode_scores_r2)
# barcode_scores_control_no_nan = filter_nans(barcode_scores_control)
# pickle.dump(barcode_scores_control_no_nan, open('./input_output_files/output/control_full_data_barcodes_no_nans.pkl', 'wb'))


barcode_scores_r1_no_nan = pickle.load(open('./input_output_files/output/r1_full_data_barcodes_no_nans.pkl', 'rb'))
barcode_scores_r2_no_nan = pickle.load(open('./input_output_files/output/r2_full_data_barcodes_no_nans.pkl', 'rb'))
barcode_scores_control_no_nan = pickle.load(open('./input_output_files/output/r2_full_data_barcodes_no_nans.pkl', 'rb'))



# ------------------------------------------------------------------------------------------------
# Make hist
# ------------------------------------------------------------------------------------------------


# print barcode_scores_r1_no_nan
r1_scores = np.array([float(tup[1]) for tup in barcode_scores_r1_no_nan.values()])
r1_scores_sorted = sorted(r1_scores, reverse = True)# counts, bin_edges = np.histogram(data, bins=num_bins, normed=True)
# print r1_scores_sorted
r1_cutoff_score = r1_scores_sorted[int(len(r1_scores_sorted) * .20)]
# print r1_cutoff_score

sns.kdeplot(r1_scores, cumulative = True)
plt.plot([r1_cutoff_score,r1_cutoff_score],[0,1],'r--', label = '20% Score Cutoff (Std Err = ' + str(r1_cutoff_score) + ')')
plt.legend(loc = 4)
plt.title('R1 Standard Error Scores (CDF)')
plt.show()


# print barcode_scores_r1_no_nan
r2_scores = np.array([float(tup[1]) for tup in barcode_scores_r2_no_nan.values()])
r2_scores_sorted = sorted(r2_scores, reverse = True)# counts, bin_edges = np.histogram(data, bins=num_bins, normed=True)
# print r1_scores_sorted
r2_cutoff_score = r2_scores_sorted[int(len(r2_scores_sorted) * .20)]
# print r1_cutoff_score

sns.kdeplot(r2_scores, cumulative = True)
plt.plot([r2_cutoff_score,r2_cutoff_score],[0,1],'r--', label = '20% Score Cutoff (Std Err = ' + str(r2_cutoff_score) + ')')
plt.legend(loc = 4)
plt.title('R2 Standard Error Scores (CDF)')
plt.show()


# # print barcode_scores_r1_no_nan
control_scores = np.array([float(tup[1]) for tup in barcode_scores_control_no_nan.values()])
control_scores_sorted = sorted(control_scores, reverse = True)# counts, bin_edges = np.histogram(data, bins=num_bins, normed=True)
# print r1_scores_sorted
control_cutoff_score = control_scores_sorted[int(len(control_scores_sorted) * .20)]
# print r1_cutoff_score

sns.kdeplot(control_scores, cumulative = True)
plt.plot([control_cutoff_score,control_cutoff_score],[0,1],'r--', label = '20% Score Cutoff (Std Err = ' + str(control_cutoff_score) + ')')
plt.legend(loc = 4)
plt.title('Control Standard Error Scores (CDF)')
plt.show()


print 'r1: ' + str(r1_cutoff_score)
print 'r2: ' + str(r2_cutoff_score)
print 'control: ' + str(control_cutoff_score)


