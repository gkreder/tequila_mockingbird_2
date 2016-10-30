import pickle
import matplotlib.pyplot as plt
import numpy as np
# import seaborn as sns
import scipy.stats

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

# R = np.corrcoef(control_fitness,weighted_avg_fitness)[0][1]
# print(scipy.stats.linregress(control_fitness, weighted_avg_fitness))[4]
slope, intercept, R, p_value, stderr = scipy.stats.linregress(control_fitness, weighted_avg_fitness)
# slope, intercept, R, p-value, stderr = scipy.stats.linregress(control_fitness, weighted_avg_fitness)
R2 = R * R

best_fit_x = np.linspace(-1,0.8,100)
best_fit_line = best_fit_x * float(R) + intercept


# --------------------------------------------------
# Scatter Plot
# --------------------------------------------------
plt.figure(figsize=[8,8])
plt.plot([-1,0.8],[-1,0.8],'r--',alpha = 1,label = "y = x")
# plt.plot(best_fit_x,best_fit_line,'r--',alpha = 0.6,label = "R2 = %.2f"%R2)
plt.scatter(control_fitness, weighted_avg_fitness, color = 'black', s = 3)
plt.xlim(-1,0.8)
plt.ylim(-1,0.8)
plt.xlabel('control_fitness',fontsize=12)
plt.ylabel('weighted_avg_fitness',fontsize = 12)
plt.legend(loc=0)
plt.savefig('control_vs_weight_avg_fitness.png',dpi=300)
plt.show()
# --------------------------------------------------



# --------------------------------------------------
# histogram
# --------------------------------------------------
# control_mean = np.mean(control_fitness)
# weighted_avg_mean = np.mean(weighted_avg_fitness)
# diffs = np.subtract(weighted_avg_fitness, control_fitness)
# diffs_mean = np.mean(diffs)
# fig, ax = plt.subplots(figsize = [8,8])
# t, prob = scipy.stats.ttest_1samp(diffs, 0.0)
# # plt.figure(figsize=[8,8])
# # plt.plot([-1,0.8],[-1,0.8],'r--',alpha = 0.6,label = "R2 = %.2f"%R2)
# # plt.hist(control_fitness, color = 'b', bins = 100, alpha = 0.8, label = 'Control')
# # plt.hist(weighted_avg_fitness, color = 'r', bins = 100, alpha = 0.5, label = 'Treatment (Weighted Average)')
# plt.hist(diffs, color = 'b', bins = 100, alpha = 0.8)
# ylim = ax.get_ylim()
# plt.plot([diffs_mean, diffs_mean], [0, ylim[1]], 'b--', alpha = 0.8, label = 'p-value (two-tailed) = ' + str(prob))
# plt.plot([0, 0], [0, ylim[1]], 'k', alpha = 0.3)
# # plt.legend()
# plt.savefig('weighted_avg_control_diff_hist.png',dpi=300)
# plt.show()
# --------------------------------------------------

