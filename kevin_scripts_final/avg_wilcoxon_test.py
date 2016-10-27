from scipy import stats
import pickle
import matplotlib.pyplot as plt

control_dict = pickle.load(open('control_filtered_AA_scores_20.pkl', 'r'))
avg_dict = pickle.load(open('weighed_avg_r1_r2_full_data.pkl', 'r'))

position_aa_dict = {}
control_position_aa_dict = {}
working_dict = avg_dict
for key in working_dict.keys():
	if key in control_dict.keys():
		position = key[0]
		aa = key[1]
		score = working_dict[key]
		control_score = control_dict[key]['score']
		if position != 0:
			if position in position_aa_dict.keys():
				position_aa_dict[position][aa] = score
				control_position_aa_dict[position][aa] = control_score
			else:
				position_aa_dict[position] = {aa : score}
				control_position_aa_dict[position] = {aa : control_score}

position_aa_list = {}
control_position_aa_list = {}

for position in position_aa_dict.keys():
	position_aa_list[position] = []
	control_position_aa_list[position] = []
	for aa in position_aa_dict[position].keys():
		exp_score = position_aa_dict[position][aa]
		control_score = control_position_aa_dict[position][aa]
		position_aa_list[position].append(exp_score)
		control_position_aa_list[position].append(control_score)

p_value_dict = {}

for position in position_aa_list.keys():
	a,p = stats.wilcoxon(position_aa_list[position], control_position_aa_list[position])
	p_value_dict[position] = p

print p_value_dict
with open("avg_weighed_p_value_dict.pkl", "wb") as f:
	pickle.dump(p_value_dict, f)
