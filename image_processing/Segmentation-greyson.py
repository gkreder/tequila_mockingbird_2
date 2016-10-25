import skimage.morphology as morph
import skimage as sk
import numpy as np
import scipy.ndimage as ndi
import skimage.measure as measure
import matplotlib.pyplot as plt
from scipy import stats
import nd2reader as nd2r
import seaborn as sns
import skimage.feature as feature
import skimage.filters as filters
import skimage.io as io
import skimage.exposure as exposure
import skimage.segmentation as seg
from os import listdir
import skimage.measure as measure
# ----------------------------------------------------------------------
from skimage.morphology import watershed, disk
from skimage.filters import sobel
from skimage import exposure
from skimage import transform as tf
from skimage.segmentation import slic, join_segmentations
from skimage import data, img_as_float
from scipy import fftpack
from skimage.filters.rank import median
from scipy import ndimage as ndi
from skimage import data
# ----------------------------------------------------------------------
import skimage.feature
import skimage.filters
import os
import sys
import skimage
import nd2reader
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.misc
import numpy as np
GIT_DIR = os.path.normpath(os.path.dirname(os.getcwd()))
PUBS_DIR = os.path.normpath(os.path.dirname(GIT_DIR))
IMAGE_DIR = os.path.normpath(PUBS_DIR + '/images_day_1/')
# ----------------------------------------------------------------------


def normalize(img):
	high = np.amax(img)
	low = np.amin(img)
	return (img - low) / (high - low)


def find_cells(img):
	strong_blur = filters.gaussian(img, 20)
	no_back = img - strong_blur
	no_back = normalize(no_back)
	equalized_no_back = exposure.equalize_hist(no_back)
	equalized_no_back = normalize(equalized_no_back)
	edges_nb = feature.canny(equalized_no_back, sigma=5)
	close_nb = ndi.binary_closing(edges_nb, structure=np.ones((3, 3)), iterations=1)
	fill_close_nb = ndi.binary_fill_holes(close_nb)
	open_fcnb = ndi.binary_opening(fill_close_nb, structure=np.ones((10, 10)))
	open_bigger = ndi.morphology.binary_dilation(open_fcnb, iterations=5)
	border = morph.binary_dilation(open_fcnb) - open_fcnb
	dist = ndi.distance_transform_edt(open_fcnb)
	local_peaks = feature.peak_local_max(dist, min_distance=12, threshold_abs=4,
											labels=open_fcnb, indices=False)
	markers = ndi.label(local_peaks)[0]
	labels = morph.watershed(-dist, markers, mask=open_bigger)
	find_boundaries = seg.find_boundaries(labels)
	return labels

IMAGE_NAME = '/Plate000_WellA12_Seq0011C5XY1.tif'
im_path = IMAGE_DIR + IMAGE_NAME
im = scipy.misc.imread(im_path)
labels = find_cells(im)
plt.imshow(im, cmap=plt.cm.gray, interpolation='nearest')
plt.contour(labels, [0.5], linewidths=1.2, colors='y')
# plt.imshow(labeled_img)
plt.show()


