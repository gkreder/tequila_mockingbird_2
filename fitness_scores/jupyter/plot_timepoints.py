import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import numpy as np

# d = pickle.load(open('./input_output_files/output/r1_full_data_barcodes.pkl', 'rb'))

time_points = np.array([0.0, 1.91, 3.328354364])
values = np.array([0.0001, 0.0004, 0.0002])
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(time_points,values)
best_fit_line = slope * time_points + intercept
# print slope * time_points


plt.plot(time_points, values, '*')
plt.plot(time_points, best_fit_line, '--')
plt.xlim(-0.5, 5)
plt.ylim(0, 0.0005)
plt.xlabel('Generations')
plt.ylabel('Fraction of Total Reads')
plt.title('Single Barcode Time Points')
# # plt.ylabels([''])
plt.show()

