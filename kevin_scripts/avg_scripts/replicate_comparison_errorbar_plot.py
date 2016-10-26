import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('r1_full_data.pkl', 'r'))
R2_dict = pickle.load(open('r2_full_data.pkl', 'r'))
# print R2_dict

fitness_rep1 = []
fitness_rep2 = []
error_rep1 = []
error_rep2 = []

for key in R1_dict[0].keys():
	if key in R2_dict[0].keys():
		fitness1 = R1_dict[0][key]['score']
		fitness2 = R2_dict[0][key]['score']
		error1 = R1_dict[0][key]['std_dev']
		error2 = R2_dict[0][key]['std_dev']
		if ~np.isnan(fitness1 * fitness2): 
			fitness_rep1.append(fitness1)
			fitness_rep2.append(fitness2)
			error_rep1.append(error1)
			error_rep2.append(error2)

R = np.corrcoef(fitness_rep1,fitness_rep2)[0][1]
R2 = R * R

plt.figure(figsize=[8,8])
plt.plot([-1,0.8],[-1,0.8],'r--',alpha = 0.6,label = "R2 = %.2f"%R2)
plt.errorbar(fitness_rep1,fitness_rep2, xerr = error_rep1, yerr = error_rep2, color = 'black', s = 1)
plt.xlim(-1,0.8)
plt.ylim(-1,0.8)
plt.xlabel('fitness_rep1',fontsize=12)
plt.ylabel('fitness_rep2',fontsize = 12)
plt.legend(loc=0)
plt.savefig('fitness_rep1_vs_rep2.png',dpi=300)

plt.show()