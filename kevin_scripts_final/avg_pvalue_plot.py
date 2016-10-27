import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('avg_weighed_p_value_dict.pkl', 'r'))

pvalues_rep1 = []

x = []

for key in R1_dict.keys():
	pvalue1 = -np.log10(R1_dict[key])
	pvalues_rep1.append(pvalue1)
	x.append(int(key))
	
cutoff = -np.log10(0.05/76)
plt.figure(figsize = [30,6])
plt.plot(x,pvalues_rep1,'k.-',markersize = 10)
plt.plot([0,80],[cutoff,cutoff],'r--')
plt.xlim(0,80)
plt.xticks(range(80),range(80),rotation = 90)
plt.xlabel('Position',fontsize=14)
plt.ylabel('-log10 (P)',fontsize=14)
plt.savefig('pvalue_positions_avg_weighed.png',dpi = 300)
# plt.errorbar(x,pvalues_avg,yerr = pvalues_std)
plt.show()