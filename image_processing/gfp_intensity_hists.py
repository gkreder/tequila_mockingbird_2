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
	filtered_intensities = filtered_intensities[int(len(filtered_intensities) * .00) : int(len(filtered_intensities) * .95)]
	return filtered_intensities

def get_expected(x):
	mu_1 = x[int(len(x) * .3)]
	mu_2 = x[int(len(x) * .7)]
	sigma_1 = np.abs(x[int(len(x) * .05)] - mu_1)
	sigma_2 = np.abs(x[int(len(x) * .95)] - mu_2)
	# sigma_1 = mu_1 / 2
	# sigma_2 = sigma_1
	# sigma_2 = abs()
	# return mu_1, sigma_1, 1, mu_2, sigma_2, 1
	# return (0.04, 0.01, 150, 0.01, 0.01, 150)
	# return (0.007, 0.003, 150, 0.017, 0.003, 150)
	return (0.007, 0.003, 150, 0.017, 0.003, 150)

def plot_intensity_hist(intensities_1, intensities_2, label_1, label_2):
	fig, ax = plt.subplots()
	# y_1, x_1 = np.histogram(intensities_1)
	weights = np.ones_like(intensities_1)/float(len(intensities_1))
	# y_1, x_1, _ = plt.hist(intensities_1, color = 'b', alpha = 0.5, bins = 100, label = label_1, weights=weights)
	y_1, x_1, _ = plt.hist(intensities_1, color = 'b', alpha = 0.5, bins = 100, label = label_1)
	# print np.sum(y_1)
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

	weights = np.ones_like(intensities_2)/float(len(intensities_2))
	# y_2, x_2, _ = plt.hist(intensities_2, color = 'r', alpha = 0.5, bins = 100, label = label_2, weights= weights)
	y_2, x_2, _ = plt.hist(intensities_2, color = 'r', alpha = 0.5, bins = 100, label = label_2)
	values_2 = y_2
	in_edges_2 = x_2
	x_2=(x_2[1:]+x_2[:-1])/2
	expected_2 = get_expected(x_2)
	params_2,cov_2=curve_fit(bimodal,x_2,y_2,expected_2)
	[mu_2_1, sigma_2_1, A_2_1, mu_2_2, sigma_2_2, A_2_2] = params_2
	sigma_2 = pylab.sqrt(pylab.diag(cov_2))
	line_2 = plt.plot(x_2,bimodal(x_2,*params_2),color='red',lw=3)
	# print '------------------------------------------'
	# print np.sum(y_1)
	# print np.amax(y_1)
	# print '------------------------------------------'

	

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
# SKG3_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA06_gfp_intensities.pkl')))
# SKG3_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB06_gfp_intensities.pkl')))
# SKG3_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC06_gfp_intensities.pkl')))
# SKG3_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD06_gfp_intensities.pkl')))
# SKG3_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE06_gfp_intensities.pkl')))
# SKG3_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF06_gfp_intensities.pkl')))
# SKG3_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG06_gfp_intensities.pkl')))
# SKG3_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH06_gfp_intensities.pkl')))

# SKG3_NO_DMSO = np.concatenate((SKG3_NO_DMSO_1, SKG3_NO_DMSO_2))
# SKG3_PERTURB_1 = np.concatenate((SKG3_PERTURB_1_1, SKG3_PERTURB_1_2))
# SKG3_PERTURB_2 = np.concatenate((SKG3_PERTURB_2_1, SKG3_PERTURB_2_2))
# SKG3_PERTURB_TOTAL = np.concatenate((SKG3_PERTURB_1, SKG3_PERTURB_2))
# SKG3_DMSO = np.concatenate((SKG3_DMSO_1, SKG3_DMSO_2))
# plot_intensity_hist(SKG3_DMSO, SKG3_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# SKG3_DMSO_EXPECTED = (0.007, 0.003, 150, 0.017, 0.003, 150)
# # SKG3_PERTURB_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# # plt.hist(SKG3_PERTURB_TOTAL, bins = 50)
# plt.show()
# ------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------
# RNR4
# ------------------------------------------------------------------------------------------------
# RNR4_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA09_gfp_intensities.pkl')))
# RNR4_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB09_gfp_intensities.pkl')))
# RNR4_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC09_gfp_intensities.pkl')))
# RNR4_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD09_gfp_intensities.pkl')))
# RNR4_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE09_gfp_intensities.pkl')))
# RNR4_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF09_gfp_intensities.pkl')))
# RNR4_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG09_gfp_intensities.pkl')))
# RNR4_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH09_gfp_intensities.pkl')))

# RNR4_NO_DMSO = np.concatenate((RNR4_NO_DMSO_1, RNR4_NO_DMSO_2))
# RNR4_PERTURB_1 = np.concatenate((RNR4_PERTURB_1_1, RNR4_PERTURB_1_2))
# RNR4_PERTURB_2 = np.concatenate((RNR4_PERTURB_2_1, RNR4_PERTURB_2_2))
# RNR4_PERTURB_TOTAL = np.concatenate((RNR4_PERTURB_1, RNR4_PERTURB_2))
# RNR4_DMSO = np.concatenate((RNR4_DMSO_1, RNR4_DMSO_2))
# plot_intensity_hist(RNR4_DMSO, RNR4_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# RNR4_DMSO_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# RNR4_PERTURB_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# # plt.hist(RNR4_PERTURB_TOTAL, bins = 50)
# # plt.show()
# ------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------
# HTB1
# ------------------------------------------------------------------------------------------------
HTB1_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA11_gfp_intensities.pkl')))
HTB1_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB11_gfp_intensities.pkl')))
HTB1_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC11_gfp_intensities.pkl')))
HTB1_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD11_gfp_intensities.pkl')))
HTB1_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE11_gfp_intensities.pkl')))
HTB1_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF11_gfp_intensities.pkl')))
HTB1_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG11_gfp_intensities.pkl')))
HTB1_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH11_gfp_intensities.pkl')))

HTB1_NO_DMSO = np.concatenate((HTB1_NO_DMSO_1, HTB1_NO_DMSO_2))
HTB1_PERTURB_1 = np.concatenate((HTB1_PERTURB_1_1, HTB1_PERTURB_1_2))
HTB1_PERTURB_2 = np.concatenate((HTB1_PERTURB_2_1, HTB1_PERTURB_2_2))
HTB1_PERTURB_TOTAL = np.concatenate((HTB1_PERTURB_1, HTB1_PERTURB_2))
HTB1_DMSO = np.concatenate((HTB1_DMSO_1, HTB1_DMSO_2))
plot_intensity_hist(HTB1_DMSO, HTB1_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
HTB1_DMSO_EXPECTED = (0.007, 0.003, 150, 0.017, 0.003, 150)
# HTB1_PERTURB_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# plt.hist(HTB1_PERTURB_TOTAL, bins = 50)
plt.show()
# ------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------
# PRE6
# ------------------------------------------------------------------------------------------------
# PRE6_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA07_gfp_intensities.pkl')))
# PRE6_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB07_gfp_intensities.pkl')))
# PRE6_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC07_gfp_intensities.pkl')))
# PRE6_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD07_gfp_intensities.pkl')))
# PRE6_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE07_gfp_intensities.pkl')))
# PRE6_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF07_gfp_intensities.pkl')))
# PRE6_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG07_gfp_intensities.pkl')))
# PRE6_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH07_gfp_intensities.pkl')))

# PRE6_NO_DMSO = np.concatenate((PRE6_NO_DMSO_1, PRE6_NO_DMSO_2))
# PRE6_PERTURB_1 = np.concatenate((PRE6_PERTURB_1_1, PRE6_PERTURB_1_2))
# PRE6_PERTURB_2 = np.concatenate((PRE6_PERTURB_2_1, PRE6_PERTURB_2_2))
# PRE6_PERTURB_TOTAL = np.concatenate((PRE6_PERTURB_1, PRE6_PERTURB_2))
# PRE6_DMSO = np.concatenate((PRE6_DMSO_1, PRE6_DMSO_2))
# plot_intensity_hist(PRE6_DMSO, PRE6_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# PRE6_DMSO_EXPECTED = (0.007, 0.003, 150, 0.017, 0.003, 150)
# # PRE6_PERTURB_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# # plt.hist(PRE6_PERTURB_TOTAL, bins = 50)
# plt.show()
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
# IHD2
# ------------------------------------------------------------------------------------------------
# IHD2_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA01_gfp_intensities.pkl')))
# IHD2_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB01_gfp_intensities.pkl')))
# IHD2_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC01_gfp_intensities.pkl')))
# IHD2_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD01_gfp_intensities.pkl')))
# IHD2_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE01_gfp_intensities.pkl')))
# IHD2_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF01_gfp_intensities.pkl')))
# IHD2_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG01_gfp_intensities.pkl')))
# IHD2_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH01_gfp_intensities.pkl')))

# IHD2_NO_DMSO = np.concatenate((IHD2_NO_DMSO_1, IHD2_NO_DMSO_2))
# IHD2_PERTURB_1 = np.concatenate((IHD2_PERTURB_1_1, IHD2_PERTURB_1_2))
# IHD2_PERTURB_2 = np.concatenate((IHD2_PERTURB_2_1, IHD2_PERTURB_2_2))
# IHD2_PERTURB_TOTAL = np.concatenate((IHD2_PERTURB_1, IHD2_PERTURB_2))
# IHD2_DMSO = np.concatenate((IHD2_DMSO_1, IHD2_DMSO_2))
# plot_intensity_hist(IHD2_DMSO, IHD2_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# # IHD2_DMSO_EXPECTED = (0.01, 0.01, 150, 0.05, 0.01, 150)
# # IHD2_PERTURB_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# # plt.hist(IHD2_PERTURB_TOTAL, bins = 50)
# plt.show()
# ------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------
# LSM1
# ------------------------------------------------------------------------------------------------
# LSM1_NO_DMSO_1 = filter_intensities(pickle.load(open('./input/WellA01_gfp_intensities.pkl')))
# LSM1_NO_DMSO_2 = filter_intensities(pickle.load(open('./input/WellB01_gfp_intensities.pkl')))
# LSM1_DMSO_1 = filter_intensities(pickle.load(open('./input/WellC01_gfp_intensities.pkl')))
# LSM1_DMSO_2 = filter_intensities(pickle.load(open('./input/WellD01_gfp_intensities.pkl')))
# LSM1_PERTURB_1_1 = filter_intensities(pickle.load(open('./input/WellE01_gfp_intensities.pkl')))
# LSM1_PERTURB_1_2 = filter_intensities(pickle.load(open('./input/WellF01_gfp_intensities.pkl')))
# LSM1_PERTURB_2_1 = filter_intensities(pickle.load(open('./input/WellG01_gfp_intensities.pkl')))
# LSM1_PERTURB_2_2 = filter_intensities(pickle.load(open('./input/WellH01_gfp_intensities.pkl')))

# LSM1_NO_DMSO = np.concatenate((LSM1_NO_DMSO_1, LSM1_NO_DMSO_2))
# LSM1_PERTURB_1 = np.concatenate((LSM1_PERTURB_1_1, LSM1_PERTURB_1_2))
# LSM1_PERTURB_2 = np.concatenate((LSM1_PERTURB_2_1, LSM1_PERTURB_2_2))
# LSM1_PERTURB_TOTAL = np.concatenate((LSM1_PERTURB_1, LSM1_PERTURB_2))
# LSM1_DMSO = np.concatenate((LSM1_DMSO_1, LSM1_DMSO_2))
# plot_intensity_hist(LSM1_DMSO, LSM1_PERTURB_TOTAL, 'DMSO CONTROL', 'PERTURB')
# # LSM1_DMSO_EXPECTED = (0.01, 0.01, 150, 0.05, 0.01, 150)
# # LSM1_PERTURB_EXPECTED = (0.025, 0.025, 150, 0.25, 0.05, 150)
# # plt.hist(LSM1_PERTURB_TOTAL, bins = 50)
# plt.show()
# ------------------------------------------------------------------------------------------------

