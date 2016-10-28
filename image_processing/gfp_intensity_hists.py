#!/usr/bin/env python2
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
#from sys import argv
from skimage import color, filters, measure, exposure, io
import skimage
from scipy import stats
import os
import pickle
import seaborn as sns
import pylab
from scipy.optimize import curve_fit


def gauss(x,mu,sigma,A):
    return A*pylab.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def filter_intensities(intensities):
	filtered_intensities = sorted(intensities)
	filtered_intensities = filtered_intensities[int(len(filtered_intensities) * .00) : int(len(filtered_intensities) * .90)]
	return filtered_intensities

def get_expected(x):
	mu_1 = x[int(len(x) * .25)]
	mu_2 = x[int(len(x) * .75)]
	sigma_1 = np.abs(x[int(len(x) * .05)] - mu_1)
	sigma_2 = np.abs(x[int(len(x) * .95)] - mu_2)
	# sigma_1 = mu_1 / 2
	# sigma_2 = sigma_1
	# sigma_2 = abs()
	return mu_1, sigma_1, 1, mu_2, sigma_2, 1

def plot_intensity_hist(intensities_1, intensities_2, label_1, label_2):
	fig, ax = plt.subplots()
	# y_1, x_1 = np.histogram(intensities_1)
	weights = np.ones_like(intensities_1)/float(len(intensities_1))
	y_1, x_1, _ = plt.hist(intensities_1, color = 'b', alpha = 0.5, bins = 100, label = label_1, weights=weights)
	print np.sum(y_1)
	values_1 = y_1
	in_edges_1 = x_1
	x_1=(x_1[1:]+x_1[:-1])/2
	expected_1 = get_expected(x_1)
	params_1,cov_1=curve_fit(bimodal,x_1,y_1,expected_1)
	[mu_1_1, sigma_1_1, A_1_1, mu_1_2, sigma_1_2, A_1_2] = params_1
	# print mu_1_1
	sigma_1 = pylab.sqrt(pylab.diag(cov_1))
	plt.plot(x_1,bimodal(x_1,*params_1),color='blue',lw=3)
	# mu_1_1_y = bimodal(x_1 * params_1)[mu_1_1]
	# print '---------------------------------'
	# print mu_1_1,mu_1_1_y, ylim_1[1]
	# print '---------------------------------'
	# print mu_1_1_y

	weights = np.ones_like(intensities_2)/float(len(intensities_1))
	y_2, x_2, _ = plt.hist(intensities_2, color = 'r', alpha = 0.5, bins = 100, label = label_2, weights= weights)
	values_2 = y_2
	in_edges_2 = x_2
	x_2=(x_2[1:]+x_2[:-1])/2
	expected_2 = get_expected(x_2)
	params_2,cov_2=curve_fit(bimodal,x_2,y_2,expected_2)
	[mu_2_1, sigma_2_1, A_2_1, mu_2_2, sigma_2_2, A_2_2] = params_2
	sigma_2 = pylab.sqrt(pylab.diag(cov_2))
	line_2 = plt.plot(x_2,bimodal(x_2,*params_2),color='red',lw=3)

	

	# mu_1_1_y = y_1[int(mu_1_1)]
	# mu_1_2_y = y_1[int(mu_1_2)]
	# mu_2_1_y = y_1[int(mu_2_1)]
	# mu_2_2_y = y_1[int(mu_2_2)]
	# ylim = ax.get_ylim()
	# # print mu_1_1_y, mu_1_2_y, mu_2_1_y, mu_2_2_y
	# # print mu_1_1, mu_1_2, mu_2_1, mu_2_2
	# # print line_2[0].get_ydata()
	# plt.plot([mu_1_1, mu_1_1], [0, ylim[1]], 'b--')
	# plt.plot([mu_1_2, mu_1_2], [0, ylim[1]], 'b--')
	# plt.plot([mu_2_1, mu_2_1], [0, ylim[1]], 'r--')
	# plt.plot([mu_2_2, mu_2_2], [0, ylim[1]], 'r--')
	# # plt.plot([mu_1_1 + sigma_1_1, mu_1_1 + sigma_1_1], [0, ylim[1]], 'b--')
	# # plt.plot([mu_1_2 + sigma_1_2, mu_1_2 + sigma_1_2], [0, ylim[1]], 'b--')
	# # plt.plot([mu_2_1 + sigma_2_1, mu_2_1 + sigma_2_1], [0, ylim[1]], 'r--')
	# # plt.plot([mu_2_2 + sigma_2_2, mu_2_2 + sigma_2_2], [0, ylim[1]], 'r--')
	


	plt.legend()
	plt.show()


# ------------------------------------------------------------------------------------------------
# SEC62
# ------------------------------------------------------------------------------------------------
# SEC62_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA05_gfp_intensities.pkl')))
# SEC62_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB05_gfp_intensities.pkl')))
# SEC62_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC05_gfp_intensities.pkl')))
# SEC62_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD05_gfp_intensities.pkl')))
# SEC62_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE05_gfp_intensities.pkl')))
# SEC62_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF05_gfp_intensities.pkl')))
# SEC62_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG05_gfp_intensities.pkl')))
# SEC62_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH05_gfp_intensities.pkl')))

# SEC62_NO_DMSO = np.concatenate((SEC62_NO_DMSO_1, SEC62_NO_DMSO_2))
# SEC62_PERTURB_1 = np.concatenate((SEC62_PERTURB_1_1, SEC62_PERTURB_1_2))
# SEC62_PERTURB_2 = np.concatenate((SEC62_PERTURB_2_1, SEC62_PERTURB_2_2))
# SEC62_PERTURB_TOTAL = np.concatenate((SEC62_PERTURB_1, SEC62_PERTURB_2))
# SEC62_DMSO = np.concatenate((SEC62_DMSO_1, SEC62_DMSO_2))
# plot_intensity_hist(SEC62_DMSO, SEC62_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# ------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------
# VMA4
# ------------------------------------------------------------------------------------------------
# VMA4_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA03_gfp_intensities.pkl')))
# VMA4_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB03_gfp_intensities.pkl')))
# VMA4_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC03_gfp_intensities.pkl')))
# VMA4_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD03_gfp_intensities.pkl')))
# VMA4_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE03_gfp_intensities.pkl')))
# VMA4_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF03_gfp_intensities.pkl')))
# VMA4_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG03_gfp_intensities.pkl')))
# VMA4_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH03_gfp_intensities.pkl')))

# VMA4_NO_DMSO = np.concatenate((VMA4_NO_DMSO_1, VMA4_NO_DMSO_2))
# VMA4_PERTURB_1 = np.concatenate((VMA4_PERTURB_1_1, VMA4_PERTURB_1_2))
# VMA4_PERTURB_2 = np.concatenate((VMA4_PERTURB_2_1, VMA4_PERTURB_2_2))
# VMA4_PERTURB_TOTAL = np.concatenate((VMA4_PERTURB_1, VMA4_PERTURB_2))
# VMA4_DMSO = np.concatenate((VMA4_DMSO_1, VMA4_DMSO_2))
# plot_intensity_hist(VMA4_DMSO, VMA4_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# SKG3
# ------------------------------------------------------------------------------------------------
SKG3_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA06_gfp_intensities.pkl')))
SKG3_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB06_gfp_intensities.pkl')))
SKG3_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC06_gfp_intensities.pkl')))
SKG3_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD06_gfp_intensities.pkl')))
SKG3_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE06_gfp_intensities.pkl')))
SKG3_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF06_gfp_intensities.pkl')))
SKG3_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG06_gfp_intensities.pkl')))
SKG3_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH06_gfp_intensities.pkl')))

SKG3_NO_DMSO = np.concatenate((SKG3_NO_DMSO_1, SKG3_NO_DMSO_2))
SKG3_PERTURB_1 = np.concatenate((SKG3_PERTURB_1_1, SKG3_PERTURB_1_2))
SKG3_PERTURB_2 = np.concatenate((SKG3_PERTURB_2_1, SKG3_PERTURB_2_2))
SKG3_PERTURB_TOTAL = np.concatenate((SKG3_PERTURB_1, SKG3_PERTURB_2))
SKG3_DMSO = np.concatenate((SKG3_DMSO_1, SKG3_DMSO_2))
plot_intensity_hist(SKG3_DMSO, SKG3_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# ------------------------------------------------------------------------------------------------



