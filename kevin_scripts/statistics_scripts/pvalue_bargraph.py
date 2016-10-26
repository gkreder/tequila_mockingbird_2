import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('R1_p_value_dict.pkl', 'r'))
R2_dict = pickle.load(open('R2_p_value_dict.pkl', 'r'))

pvalues_rep1 = []
pvalues_rep2 = []
pvalues_avg = []
pvalues_std = []
x = []
for key in R1_dict.keys():
	if key in R2_dict.keys():
		pvalue1 = -np.log10(R1_dict[key])
		pvalue2 = -np.log10(R2_dict[key])
		pvalues_rep1.append(pvalue1)
		pvalues_rep2.append(pvalue2)
		average = (pvalue1 + pvalue2) / 2.0
		std = np.std([pvalue1, pvalue2])
		pvalues_avg.append(average)
		pvalues_std.append(std)
		x.append(int(key))
cutoff = -np.log10(0.05/76)
plt.figure(figsize = [30,6])
plt.plot(x,pvalues_rep1,'k.-',markersize = 10)
plt.plot([0,80],[cutoff,cutoff],'r--')
plt.xlim(0,80)
plt.xticks(range(80),range(80),rotation = 90)
plt.xlabel('Position',fontsize=14)
plt.ylabel('-log10 (P)',fontsize=14)
plt.savefig('pvalue_positions_rep1.png',dpi = 500)
# plt.errorbar(x,pvalues_avg,yerr = pvalues_std)
plt.show()