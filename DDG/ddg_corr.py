import csv
import pickle
import scipy.stats



ddg_affinities_OTU = {}
with open('uby_1ubq_monomer.tsv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row_number, row in enumerate(spamreader):
		# print row
		# if row_number != 0:
			# print row[1].split(' ')[1][-1]
		if row_number == 0:
			continue
		AA = row[1].split(' ')[1][-1]
		pos = int(row[1].split(' ')[1][1 : -1])
		# print AA, pos
		index_tup = (pos, AA)
		# # print index_tup
		# # print row[3]
		# # (For monomer)
		if row[2] != 'None':
		# # if row[3] != 'None': (for non monomer)
			ddg_affinities_OTU[index_tup] = float(row[2])


r1_differences_scores = pickle.load(open('../fitness_scores/jupyter/input_output_files/output/r1_differences_scores.pkl', 'rb'))
r2_differences_scores = pickle.load(open('../fitness_scores/jupyter/input_output_files/output/r2_differences_scores.pkl', 'rb'))

common_keys_r1 = set(ddg_affinities_OTU.keys()).intersection(set(r1_differences_scores.keys()))
common_keys_r2 = set(ddg_affinities_OTU.keys()).intersection(set(r2_differences_scores.keys()))


r1_corr_scores_temp = {}
for pos, AA in common_keys_r1:
	key = (pos, AA)
	if pos in r1_corr_scores_temp:
		r1_corr_scores_temp[pos]['ddgs'].append(ddg_affinities_OTU[key])
		r1_corr_scores_temp[pos]['scores'].append(r1_differences_scores[key])
	else:
		r1_corr_scores_temp[pos] = {'ddgs' : [ddg_affinities_OTU[key]], 'scores' : [r1_differences_scores[key]]}

r2_corr_scores_temp = {}
for pos, AA in common_keys_r2:
	key = (pos, AA)
	if pos in r2_corr_scores_temp:
		r2_corr_scores_temp[pos]['ddgs'].append(ddg_affinities_OTU[key])
		r2_corr_scores_temp[pos]['scores'].append(r2_differences_scores[key])
	else:
		r2_corr_scores_temp[pos] = {'ddgs' : [ddg_affinities_OTU[key]], 'scores' : [r2_differences_scores[key]]}

r1_corr_scores = {}
for pos in r1_corr_scores_temp.keys():
	# key = (pos, AA)
	temp_x = r1_corr_scores_temp[pos]['scores']
	temp_y = r1_corr_scores_temp[pos]['ddgs']
	slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(temp_x, temp_y)
	r1_corr_scores[pos] = rvalue
pickle.dump(r1_corr_scores, open('./output/r1_ddg_affinity_mono_corr.pkl', 'wb'))

r2_corr_scores = {}
for pos in r2_corr_scores_temp.keys():
	# key = (pos, AA)
	temp_x = r2_corr_scores_temp[pos]['scores']
	temp_y = r2_corr_scores_temp[pos]['ddgs']
	slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(temp_x, temp_y)
	r2_corr_scores[pos] = rvalue
pickle.dump(r2_corr_scores, open('./output/r2_ddg_affinity_mono_corr.pkl', 'wb'))

avg_corr_scores = {}
for pos in r1_corr_scores.keys():
	avg_corr = (r1_corr_scores[pos] + r2_corr_scores[pos]) / 2.0
	avg_corr_scores[pos] = avg_corr
pickle.dump(avg_corr_scores, open('./output/avg_ddg_affinity_mono_corr.pkl', 'wb'))


# print r2_corr_scores
	# print pos




