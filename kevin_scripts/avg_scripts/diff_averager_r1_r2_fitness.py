import pickle
import numpy as np

R1_dict = pickle.load(open('r1_differences_scores.pkl', 'r'))
R2_dict = pickle.load(open('r2_differences_scores.pkl', 'r'))

avg_r1_r2_full_data = {}

for key in R1_dict.keys():
	if key in R2_dict.keys():
		fitness1 = R1_dict[key]
		fitness2 = R2_dict[key]
		if ~np.isnan(fitness1*fitness2):
			avg_r1_r2_full_data[key] = (fitness1 + fitness2) / 2.0
		elif ~np.isnan(fitness1) and np.isnan(fitness2):
			avg_r1_r2_full_data[key] = fitness1
		elif np.isnan(fitness1) and ~np.isnan(fitness2):
			avg_r1_r2_full_data[key] = fitness2

print(avg_r1_r2_full_data)

with open('diff_avg_r1_r2_full_data.pkl', 'w') as f:
	pickle.dump(avg_r1_r2_full_data, f)