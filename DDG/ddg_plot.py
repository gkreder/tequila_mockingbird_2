import pickle
import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns

# R1_dict = pickle.load(open('avg_weighed_p_value_dict.pkl', 'r'))
ub_seq = 'QIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
ub_seq_list = list(ub_seq)
corr_dict = pickle.load(open('./output/avg_ddg_affinity_mono_corr.pkl', 'rb'))


# print corr_dict

for index, c in enumerate(ub_seq_list):
    ub_seq_list[index] = str(index + 2) + ' (' + c + ')'
    # ub_seq_list[index] = 'blah blah blah'
loc_labels = ub_seq_list


corr_scores = []

positions = []

for key in corr_dict.keys():
	# if key > 1 and key < 77:
	corr_score = corr_dict[key]
	corr_scores.append(corr_score)
	positions.append(int(key))

# # print x
# # print corr_scores
# 34, 35, 40, 27

corr_scores_mean = np.mean(corr_scores)
corr_scores_std = np.std(corr_scores)
	
# cutoff = -np.log10(0.05/76)
plt.figure(figsize = [30,6])
plt.plot(positions,corr_scores,'k.',markersize = 10)
plt.plot([0, 77], [0, 0], 'k--')
plt.plot([0, 77], [corr_scores_mean, corr_scores_mean])
# plt.plot([0, 77], [2.0 * corr_scores_std, 2.0 * corr_scores_std], 'k--')
# plt.plot([0, 77], [-2.0 * corr_scores_std, -2.0 * corr_scores_std], 'k--')

# for pos in positions:
	# plt.plot(pos, corr_scores[pos - 2], 'r')

# plt.plot([2,76],[cutoff,cutoff],'r--')
plt.xlim(2,76)
# # plt.xticks(range(2,77),range(2,77),rotation = 90)
plt.xticks(range(2,77), loc_labels, rotation = 90)

plt.xlabel('Position',fontsize=14)
plt.ylabel('r',fontsize=14)
# plt.savefig('ddg_scores_corr.png',dpi = 300)
# # plt.errorbar(x,pvalues_avg,yerr = pvalues_std)


# plt.hist(corr_scores, bins = 25)
# plt.()




plt.show()