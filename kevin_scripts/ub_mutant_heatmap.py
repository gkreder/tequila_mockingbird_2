import pickle
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import sys

arg,arg1 = sys.argv
FILE = arg1
fitness_scores = pickle.load(open(FILE, 'r'))

# test_data = [{'A' : (0.1 , 0.1, 0.2, 0.15), 'V' : (0.4, 0.4, 0.5, 0.45)}, {'A' : (0.4, 0.2, 0.3, 0.25), 'V' : (0.4, 0.5, 0.6, 0.55)}]

averagefitnessdata_list = []

for position in fitness_scores:
	averageposition_dict = {}
	for aa in position.keys():
		value = position[aa]
		# averageposition_dict[aa] = value[3] - value[0]
		averageposition_dict[aa] = value[3]
	averagefitnessdata_list.append(averageposition_dict)
plt.figure(figsize=[20,8])


sns.heatmap(pd.DataFrame(averagefitnessdata_list).T, cmap = plt.cm.Reds)
plt.show()