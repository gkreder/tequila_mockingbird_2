import pickle
import numpy as np

R1_dict = pickle.load(open('r1_filtered_AA_scores_20.pkl', 'r'))
R2_dict = pickle.load(open('r2_filtered_AA_scores_20.pkl', 'r'))

avg_r1_r2_full_data = {}

def weighed_avg(fit1, fit2, std_dev1, std_dev2):
	weight1 = std_dev1 / (std_dev1 + std_dev2)
	weight2 = std_dev2 / (std_dev1 + std_dev2)
	avg = fit1 * weight1 + fit2 * weight2
	return avg

for key in R1_dict.keys():
	if key in R2_dict.keys():
		fitness1 = R1_dict[key]['score']
		fitness2 = R2_dict[key]['score']
		std_dev1 = R1_dict[key]['std_dev']
		std_dev2 = R2_dict[key]['std_dev']
		if ~np.isnan(fitness1*fitness2):
			avg_r1_r2_full_data[key] = weighed_avg(fitness1, fitness2, std_dev1, std_dev2)
		elif ~np.isnan(fitness1) and np.isnan(fitness2):
			avg_r1_r2_full_data[key] = fitness1
		elif np.isnan(fitness1) and ~np.isnan(fitness2):
			avg_r1_r2_full_data[key] = fitness2

print(avg_r1_r2_full_data)

with open('weighed_avg_r1_r2_full_data.pkl', 'w') as f:
	pickle.dump(avg_r1_r2_full_data, f)