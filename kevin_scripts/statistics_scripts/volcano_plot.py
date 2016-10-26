import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('R1_p_value_dict.pkl', 'r'))
R2_dict = pickle.load(open('R2_p_value_dict.pkl', 'r'))
R1_difference_dict = pickle.load(open('r1_differences_scores.pkl', 'r'))
R2_difference_dict = pickle.load(open('r2_differences_scores.pkl', 'r'))

position_difference_dict = {}

for key in R1_difference_dict.keys():
	position = key[0]
	fitness_difference = R1_difference_dict[key]
	if position in position_difference_dict.keys():
		position_difference_dict[position].append(fitness_difference)
	else:
		position_difference_dict[position] = [fitness_difference]

position_difference_avg_dict = {}

for position in position_difference_dict.keys():
	position_difference_avg_dict[position] = np.mean(position_difference_dict[position])

x_differences = []
y_pvalues = []

for position in position_difference_avg_dict.keys():
	if position != 0:
		x_differences.append(position_difference_avg_dict[position])
		y_pvalues.append(-np.log10(R1_dict[position]))
p_cutoff = -np.log10(0.05/76)
fitness_diff_cutoff = -0.14

plt.figure(figsize=[10,6])
plt.plot([-0.3,0.1],[p_cutoff,p_cutoff],'r--',label = 'P cutoff')
plt.plot([fitness_diff_cutoff,fitness_diff_cutoff],[0,4.5],'b--',label = 'fitness difference cutoff')
plt.xlim(-0.3,0.1)
plt.ylim(0,4.5)
plt.scatter(x_differences, y_pvalues, color = 'black')
plt.legend(loc=0)
plt.savefig('volcanoplot_rep1.png',dpi = 500)
plt.show()