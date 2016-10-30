import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats


DMSO_eccentricities = pickle.load(open('DMSO_eccentricities.pkl', 'rb'))
PERTURB_eccentriticies = pickle.load(open('PERTURB_eccentricities.pkl', 'rb'))

print DMSO_eccentricities



# data = np.loadtxt('Filename.txt')

# Choose how many bins you want here
num_bins = 50

# Use the histogram function to bin the data
DMSO_counts, DMSO_bin_edges = np.histogram(DMSO_eccentricities, bins=num_bins, normed=True)
PERTURB_counts, PERTURB_bin_edges = np.histogram(PERTURB_eccentriticies, bins=num_bins, normed=True)


DMSO_counts=np.array(DMSO_counts) / float(len(DMSO_counts) + 1)
PERTURB_counts=np.array(PERTURB_counts) / float(len(PERTURB_counts) + 1)

# Now find the cdf
DMSO_cdf = np.cumsum(DMSO_counts)
PERTURB_cdf = np.cumsum(PERTURB_counts)

max_dist = 0.0
max_index = 0
for index,val in enumerate(DMSO_cdf):
	if abs(PERTURB_cdf[index] - val) > max_dist:
		max_index = index
		max_dist = abs(PERTURB_cdf[index] - val)


# ks test
KS_stat, pval = scipy.stats.ks_2samp(DMSO_cdf, PERTURB_cdf)
print KS_stat, pval


# And finally plot the cdf
plt.plot(DMSO_bin_edges[1:], DMSO_cdf, label = 'DMSO CONTROL', color = 'b')
plt.plot(PERTURB_bin_edges[1:], PERTURB_cdf, label = 'PERTURB', color = 'r')
plt.plot([PERTURB_bin_edges[max_index + 1], PERTURB_bin_edges[max_index + 1]], 
	[DMSO_cdf[max_index], PERTURB_cdf[max_index]], 'k--' , label = 'KS Distance')
plt.legend(loc=2)
plt.title('Cell Eccentricity CDF')
plt.xlabel('Eccentricity')
plt.ylim(0,1)
plt.show()