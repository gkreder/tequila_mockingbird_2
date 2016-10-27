import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('r1_filtered_AA_scores_20.pkl', 'r'))#[0]
R2_dict = pickle.load(open('r2_filtered_AA_scores_20.pkl', 'r'))#[0]

fitness_rep1 = []
fitness_rep2 = []
error_rep1 = []
error_rep2 = []

for key in R1_dict.keys():
	if key in R2_dict.keys():
		fitness1 = R1_dict[key]['score']
		fitness2 = R2_dict[key]['score']
		error1 = R1_dict[key]['std_dev']
		error2 = R2_dict[key]['std_dev']
		if ~np.isnan(fitness1 * fitness2): 
			fitness_rep1.append(fitness1)
			fitness_rep2.append(fitness2)
			error_rep1.append(error1/10.0)
			error_rep2.append(error2/10.0)

R = np.corrcoef(fitness_rep1,fitness_rep2)[0][1]
R2 = R * R

plt.figure(figsize=[8,8])
plt.plot([-1,0.8],[-1,0.8],'r--',alpha = 0.6,label = "R2 = %.2f"%R2)
plt.errorbar(fitness_rep1,fitness_rep2, xerr = error_rep1, yerr = error_rep2,linestyle="None",color='black', fmt='.', markersize='4', ecolor='grey',capsize=1, elinewidth=1)
plt.xlim(-1,0.8)
plt.ylim(-1,0.8)
plt.xlabel('fitness_rep1',fontsize=12)
plt.ylabel('fitness_rep2',fontsize = 12)
plt.legend(loc=0)
plt.savefig('fitness_rep1_vs_rep2_errorbar.png',dpi=300)

plt.show()