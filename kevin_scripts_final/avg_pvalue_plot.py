import pickle
import matplotlib.pyplot as plt
import numpy as np

R1_dict = pickle.load(open('avg_weighed_p_value_dict.pkl', 'r'))
ub_seq = 'QIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
ub_seq_list = list(ub_seq)
for index, c in enumerate(ub_seq_list):
    ub_seq_list[index] = str(index + 2) + ' (' + c + ')'
    # ub_seq_list[index] = 'blah blah blah'
loc_labels = ub_seq_list


pvalues_rep1 = []

x = []

for key in R1_dict.keys():
	# if key > 1 and key < 77:
	pvalue1 = -np.log10(R1_dict[key])
	pvalues_rep1.append(pvalue1)
	x.append(int(key))

# print x
# print pvalues_rep1
	
cutoff = -np.log10(0.05/76)
plt.figure(figsize = [30,6])
plt.plot(x,pvalues_rep1,'k.-',markersize = 10)
# plt.plot([0,80],[cutoff,cutoff],'r--')
# plt.xlim(0,80)
# plt.xticks(range(80),range(80),rotation = 90)

plt.plot([2,76],[cutoff,cutoff],'r--')
plt.xlim(2,76)
# plt.xticks(range(2,77),range(2,77),rotation = 90)
plt.xticks(range(2,77), loc_labels, rotation = 90)

plt.xlabel('Position',fontsize=14)
plt.ylabel('-log10 (P)',fontsize=14)
plt.savefig('pvalue_positions_avg_weighed.png',dpi = 300)
# plt.errorbar(x,pvalues_avg,yerr = pvalues_std)
plt.show()