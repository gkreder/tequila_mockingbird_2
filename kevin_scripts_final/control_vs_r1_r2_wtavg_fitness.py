import pickle
import matplotlib.pyplot as plt
import numpy as np

control_fitness_dict = pickle.load(open('control_filtered_AA_scores_20.pkl', 'r'))#[0]
weighted_avg_fitness_dict = pickle.load(open('weighed_avg_r1_r2_full_data.pkl', 'r'))#[0]

control_fitness = []
weighted_avg_fitness = []

for key in control_fitness_dict.keys():
	if key in weighted_avg_fitness_dict.keys():
		fitness1 = control_fitness_dict[key]['score']
		fitness2 = weighted_avg_fitness_dict[key]
		if ~np.isnan(fitness1 * fitness2): 
			control_fitness.append(fitness1)
			weighted_avg_fitness.append(fitness2)

R = np.corrcoef(control_fitness,weighted_avg_fitness)[0][1]
R2 = R * R

plt.figure(figsize=[8,8])
plt.plot([-1,0.8],[-1,0.8],'r--',alpha = 0.6,label = "R2 = %.2f"%R2)
plt.scatter(control_fitness, weighted_avg_fitness)
plt.xlim(-1,0.8)
plt.ylim(-1,0.8)
plt.xlabel('control_fitness',fontsize=12)
plt.ylabel('weighted_avg_fitness',fontsize = 12)
plt.legend(loc=0)
plt.savefig('control_vs_weight_avg_fitness.png',dpi=300)

plt.show()