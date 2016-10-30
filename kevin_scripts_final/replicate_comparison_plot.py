import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

R1_dict = pickle.load(open('r1_full_data.pkl', 'r'))
R2_dict = pickle.load(open('r2_full_data.pkl', 'r'))
# print R2_dict

fitness_rep1 = []
fitness_rep2 = []

for key in R1_dict[0].keys():
	if key in R2_dict[0].keys():
		fitness1 = R1_dict[0][key]['score']
		fitness2 = R2_dict[0][key]['score']
		if ~np.isnan(fitness1 * fitness2): 
			fitness_rep1.append(fitness1)
			fitness_rep2.append(fitness2)


slope, intercept, rval, pval, stderr = scipy.stats.linregress(fitness_rep1, fitness_rep2)
x = np.linspace(-1,0.8, 100)
best_fit_line = float(slope) * x + intercept

R = np.corrcoef(fitness_rep1,fitness_rep2)[0][1]
R2 = R * R
plt.figure(figsize=[8,8])
plt.scatter(fitness_rep1,fitness_rep2,color = 'black',s = 3)

# plt.plot([-1,0.8],[-1,0.8],'b--',alpha = 0.6,label = "y = x")
plt.plot(x, best_fit_line,'r',alpha = 0.6,label = "R2 = %.2f"%R2)

plt.xlim(-1,0.8)
plt.ylim(-1,0.8)
plt.xlabel('fitness_rep1',fontsize=12)
plt.ylabel('fitness_rep2',fontsize = 12)
plt.legend(loc=0)
plt.savefig('fitness_rep1_vs_rep2.png',dpi=300)

plt.show()