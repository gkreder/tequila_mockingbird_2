import pickle

control_dict = pickle.load(open('control_filtered_AA_scores_20.pkl', 'r'))
weighed_avg_dict = pickle.load(open('weighed_avg_r1_r2_full_data.pkl', 'r'))

weighed_avg_fitness_difference_dict = {}

for key in control_dict.keys():
	if key in weighed_avg_dict.keys():
		control_fitness = control_dict[key]['score']
		weighed_avg_fitness = weighed_avg_dict[key]
		weighed_avg_fitness_difference_dict[key] = weighed_avg_fitness - control_fitness

print weighed_avg_fitness_difference_dict

with open('weighed_avg_fitness_difference_dict.pkl', 'w') as f:
	pickle.dump(weighed_avg_fitness_difference_dict, f)