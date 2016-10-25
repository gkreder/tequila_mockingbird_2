import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

dmso_intensities = pickle.load(open('dmso_intensities.pkl'))
perturb_intensities = pickle.load(open('perturb_intensities.pkl'))

# # print dmso_intensities
# plt.hist(dmso_intensities, 'r')
# plt.hist(perturb_intensities, 'b',  alpha = 0.5)
# plt.show()

plt.figure()
plt.title('RNR4 Integrated GFP Intensity Per Cell')
sns.kdeplot(np.array(dmso_intensities),shade=True,color='r', label = 'dmso control')
sns.kdeplot(np.array(perturb_intensities),shade=True,color='b', label = 'ampho B perturb')
plt.legend()
plt.show()