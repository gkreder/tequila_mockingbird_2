import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('avg_p_value_dict.pkl', 'r'))
R1_difference_dict = pickle.load(open('diff_avg_r1_r2_full_data.pkl', 'r'))

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
positions = []

for position in position_difference_avg_dict.keys():
	if position != 0:
		positions.append(position)
		x_differences.append(position_difference_avg_dict[position])
		y_pvalues.append(-np.log10(R1_dict[position]))
p_cutoff = -np.log10(0.05/76)
fitness_diff_cutoff = -0.15

plt.figure(figsize=[10,8])
plt.plot([-0.3,0.1],[p_cutoff,p_cutoff],'r--',label = 'P cutoff')
plt.plot([fitness_diff_cutoff,fitness_diff_cutoff],[0,4.5],'b--',label = 'fitness difference cutoff')
plt.xlim(-0.3,0.1)
plt.ylim(0,4.5)
plt.xlabel('delta_fitness')
plt.ylabel('-log10(P)')
plt.scatter(x_differences, y_pvalues, color = 'grey')
plt.legend(loc=0)
for i in range(len(positions)):
	if (x_differences[i] <= fitness_diff_cutoff) and (y_pvalues[i] >= p_cutoff):
		plt.scatter(x_differences[i],y_pvalues[i],color = 'black')
		plt.annotate(positions[i], (x_differences[i], y_pvalues[i]),xytext=(0,6), textcoords='offset points')
		if positions[i] in [34,35,37,40]:
				plt.scatter(x_differences[i],y_pvalues[i],color = 'red')

plt.savefig('volcanoplot_avg.png',dpi = 500)
plt.show()