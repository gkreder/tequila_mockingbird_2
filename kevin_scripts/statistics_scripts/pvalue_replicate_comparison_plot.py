import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('R1_p_value_dict.pkl', 'r'))
R2_dict = pickle.load(open('R2_p_value_dict.pkl', 'r'))
# print R2_dict

pvalues_rep1 = []
pvalues_rep2 = []

for key in R1_dict.keys():
	if key in R2_dict.keys():
		pvalue1 = -np.log10(R1_dict[key])
		pvalue2 = -np.log10(R2_dict[key])
		pvalues_rep1.append(pvalue1)
		pvalues_rep2.append(pvalue2)

R = np.corrcoef(pvalues_rep1, pvalues_rep2)[0][1]
R2 = R * R

plt.figure(figsize=[8,8])
plt.scatter(pvalues_rep1,pvalues_rep2, color = 'black',s = 10)
plt.plot([-1,5],[-1,5],'r--',alpha = 0.6,label = 'R2 = %.2f'%R2)
plt.xlim(-1,5)
plt.ylim(-1,5)
plt.xlabel('-log10(P)_rep1',fontsize=12)
plt.ylabel('-log10(P)_rep2',fontsize = 12)
plt.savefig('pvalue_rep1_vs_rep2.png',dpi=300)
plt.legend(loc=0)
plt.show()