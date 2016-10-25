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
import seaborn as sns
import pickle
GIT_DIR = os.path.normpath(os.path.dirname(os.getcwd()))
PUBS_DIR = os.path.normpath(os.path.dirname(GIT_DIR))
global IMAGE_DIR
IMAGE_DIR = os.path.normpath(PUBS_DIR + '/images_day_1/')
# sys.path.append(IMAGE_DIR)
# print 'SYSTEM PATH'
# print sys.path
# ----------------------------------------------------------------------
def segment_BF(img):

	sobel_edges = sobel(img)
	global_thresh = skimage.filters.threshold_otsu(img)
	binary_global = img > global_thresh
	elevation_map = binary_global

	# img_norm = skimage.img_as_float(img)
	img_norm = img / img.astype(np.uint16)
	# img_norm = img / img.astype(np.float64)
	markers = np.zeros_like(img)
	markers[img > np.mean(img) - 200] = 1
	markers[img < np.mean(img) - 200] = 2
	# markers = skimage.morphology.closing(markers)
	# markers = ndi.binary_fill_holes(markers - 1)

	# ------------------------------------------------------------------------
	# We should try changing the marker matrix eventually
	# ------------------------------------------------------------------------
	# local_peaks = skimage.feature.peak_local_max(img_norm, min_distance = 10)
	# # , threshold_abs = 0.3
	# for peak_row, peak_column in local_peaks:
	# 	markers[peak_row, peak_column] = 1
	# # markers = markers + 1
	# ------------------------------------------------------------------------

	segmentation = watershed(elevation_map, markers)

	segmentation = skimage.morphology.closing(segmentation)	

	segmentation = ndi.binary_fill_holes(segmentation - 1)
	segmentation_first_filter = skimage.morphology.remove_small_objects(segmentation, 200)
	segmentation_second_filter = skimage.morphology.remove_small_objects(segmentation, 800)
	segmentation_second_filter = np.invert(segmentation_second_filter)
	segmentation_final = np.multiply(segmentation_first_filter, segmentation_second_filter)
	labeled_img, num_objects = ndi.label(segmentation_final)
	image_label_overlay = skimage.color.label2rgb(labeled_img, image=img)

	# plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
	# plt.contour(segmentation_final, [0.5], linewidths=1.2, colors='y')

	# plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
	# plt.contour(markers, [0.5], linewidths=1.2, colors='y')

	# plt.imshow(image_label_overlay)
	# plt.imshow(markers, cmap = plt.cm.gray)
	# plt.show()

	return labeled_img, num_objects
	# return labeled_img[0], num_objects


def segment_DAPI(img):
	marker_mat = np.zeros_like(img)
	global_thresh = skimage.filters.threshold_otsu(img)
	binary_global = img > global_thresh
	sobel_edges = sobel(img)

	markers = binary_global + 1
	segmentation = watershed(sobel_edges, markers)
	segmentation = ndi.binary_fill_holes(segmentation - 1)

	labeled_img, num_objects = ndi.label(segmentation)
	image_label_overlay = skimage.color.label2rgb(labeled_img, image=img)

	# print object_labels


	# plt.imshow(segmentation, cmap = plt.cm.gray) 
	# plt.imshow(img, cmap = plt.cm.gray)
	# plt.contour(segmentation, [0.5], linewidths=1.2, colors='y')
	# plt.show()

	return labeled_img, num_objects

def find_nuclei_xy(img, dapi_cells):

	# ------------------------------------------------------------
	# Total Image Maxima
	# ------------------------------------------------------------
	# img = skimage.img_as_float(img)
	img = img / np.amax(img).astype(np.uint16)

	# local_peaks = skimage.feature.peak_local_max(img, min_distance = 10, threshold_abs = 0.1)
	local_peaks = skimage.feature.peak_local_max(img, min_distance = 10)
	# print local_peaks
	nuclei_image = np.zeros_like(img)
	nuclei = []
	for row, column in local_peaks:
		nuclei_image[row, column] = 1
		nuclei.append([row, column])

	nuclei_image = nuclei_image * dapi_cells

	# plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
	# plt.contour(nuclei_image, [0.5], linewidths=1.2, colors='y')
	# plt.contour(dapi_cells, [0.5], linewidths=1.2, colors='y')
	# plt.imshow(nuclei_image, cmap = plt.cm.gray)
	# plt.show()
	# print nuclei
	# ------------------------------------------------------------


	# ----------------------------------------------------------------
	# Single cell masking
	# ----------------------------------------------------------------
	# object_labels = get_object_labels(dapi_cells)	
	# nuclei_labeled = []
	# nuclei = []
	# for label in object_labels:
	# # label = 500
	# 	masked_image = np.zeros_like(dapi_cells)
	# 	img_copy = np.copy(img)
	# 	img_copy[dapi_cells != label] = 0
	# 	local_peaks = skimage.feature.peak_local_max(img_copy, num_peaks = 1)
	# 	if len(local_peaks) > 0:
	# 		peak_row, peak_column = local_peaks[0]
	# 		nuclei_labeled.append(([peak_row, peak_column], label))
	# 		nucleu.append([peak_row, peak_column])
	# 		# img[peak_row, peak_column] = 0.0

	# # plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
	# # plt.show()
	# ----------------------------------------------------------------

	# masked_image[dapi_cells == label] = 1

		

	# print dapi_cells




	# image_label_overlay = skimage.color.label2rgb(dapi_cells, image=img)
	# plt.imshow(image_label_overlay, cmap = plt.cm.gray)

	# plt.imshow(img, cmap=plt.cm.gray, interpolation='nearest')
	# plt.contour(segmentation_final, [0.5], linewidths=1.2, colors='y')


	# plt.imshow(img_copy, cmap = plt.cm.gray)
	return nuclei

def get_masked_cell_image(intensity_image, label_image, label_index):
	mask = np.zeros_like(label_image)
	mask[label_image == label_index] = 1
	masked_image = mask * intensity_image
	return masked_image



def integrated_intensity(intensity_img, labeled_cells):
	labeled_cells_binary = labeled_cells > 0

	intensity_img_masked = labeled_cells_binary * intensity_img
	intensity_img_masked_norm = skimage.img_as_float(intensity_img_masked)
	# intensity_img_masked_norm = intensity_img_masked / intensity_img_masked.astype(np.uint16)
	object_labels = get_object_labels(labeled_cells)

	# plt.imshow(intensity_img_masked)
	# plt.show()
	regionprops = skimage.measure.regionprops(labeled_cells, intensity_image = intensity_img_masked_norm)
	cell_intensities = []
	for region in regionprops:
		masked_image = get_masked_cell_image(intensity_img_masked_norm, labeled_cells, region.label)
		cell_intensity_sum = np.sum(masked_image)
		cell_intensity_normalized = cell_intensity_sum / region.area
		cell_intensities.append(cell_intensity_normalized)
		# print cell_intensity_normalized
	# plt.imshow(masked_image, cmap = plt.cm.gray)
	# plt.show()

	# plt.hist(cell_intensities, bins = 30)
	# plt.show()
	return cell_intensities
		

def find_nuclei_boundaries(labeled_img, img):
	masked_cells_img = labeled_img * img
	auto = skimage.filters.rank.autolevel(masked_cells_img, disk(5))

	plt.contour(auto, [0.5], linewidths=1.2, colors='y')
	plt.imshow(img, cmap = plt.cm.gray)
	plt.show()
	#def get_masked_cell_image(intensity_image, label_image, label_index):



def get_object_labels(labeled_img):
	object_labels = np.unique(labeled_img)
	object_labels = object_labels[1 : ]
	return object_labels

def get_nd2_dict(img_path):

	print img_path in os.listdir(IMAGE_DIR)
	# print os.listdir(IMAGE_DIR)

	im_name = img_path.replace('.nd2', '')
	im_name_pickle = im_name + '.pkl'

	if im_name_pickle not in os.listdir(IMAGE_DIR):

		print '---------------------------------'
		print 'CREATING NEW PICKLE FILE'
		print '---------------------------------'

		img_dict = {'DAPI' : {},
					'FITC' : {},
					'FITClong' : {},
					'CY3' : {},
					'BF-Cy3' : {}}

		nd2 = nd2reader.Nd2(IMAGE_DIR+'/'+img_path)
		# nd2 = nd2reader.Nd2(img_path)

		print nd2.fields_of_view


		for i in range(len(nd2)):
			channel = nd2[i].channel
			fov = nd2[i].field_of_view
			# img_dict[channel][fov] = nd2[i]
			img_dict[channel][fov] = nd2[i].astype(np.uint16) 
			# / np.amax(nd2[i]).astype(np.uint16)

		pickle.dump(img_dict, open(IMAGE_DIR + '/' + im_name_pickle, 'w'))
	else:
		print '---------------------------------'
		print 'LOADING FROM PICKLE FILE'
		print '---------------------------------'
		img_dict = pickle.load( open(IMAGE_DIR + '/' + im_name_pickle, 'r' ))
		print '...done'

	return img_dict


	# for channel in img_dict:
		# for image in nd2.select(Channel = channel):
			# print image

# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------

# # ------------------------------------------------
# # DAPI Segmentation
# # ------------------------------------------------
# # IMAGE_NAME = '/Plate000_WellA12_Seq0011C1XY1.tif'
# # im_path = IMAGE_DIR + IMAGE_NAME
# # im = scipy.misc.imread(im_path)
# # # print im
# # labeled_img, num_objects = segment_DAPI(im)
# # # plt.imshow(labeled_img)
# # # plt.show()
# # ------------------------------------------------
# # ------------------------------------------------


# # ------------------------------------------------
# # BF
# # ------------------------------------------------
# # IMAGE_NAME = '/Plate000_WellA12_Seq0011C5XY1.tif'
# # im_path = IMAGE_DIR + IMAGE_NAME
# # im = scipy.misc.imread(im_path)
# # # print im
# # segment_BF(im)
# # # plt.imshow(im)
# # # plt.show()
# # ------------------------------------------------
# # ------------------------------------------------

# # ------------------------------------------------
# # Find nuclei
# # ------------------------------------------------
# # IMAGE_NAME = '/Plate000_WellA12_Seq0011C1XY1.tif'
# # im_path = IMAGE_DIR + IMAGE_NAME
# # im = scipy.misc.imread(im_path)
# # labeled_img, num_objects = segment_DAPI(im)
# # # print im
# # find_nuclei_xy(im, labeled_img)
# # # # plt.imshow(labeled_img)
# # # # plt.show()
# # ------------------------------------------------
# # ------------------------------------------------



# # ------------------------------------------------
# # Integrated GFP intensity
# # ------------------------------------------------
IMAGE_NAME_PERTURB_BF = '/Plate000_WellE09_Seq0056XY1C5.tif'
IMAGE_NAME_DMSO_BF = '/Plate000_WellC09_Seq0032XY1C5.tif'
im_path_perturb_bf = IMAGE_DIR + IMAGE_NAME_PERTURB_BF
im_path_dmso_bf = IMAGE_DIR + IMAGE_NAME_DMSO_BF
im_perturb_bf = scipy.misc.imread(im_path_perturb_bf)
im_dmso_bf = scipy.misc.imread(im_path_dmso_bf)

# plt.imshow(im_perturb_bf)
# plt.show()


labeled_img_perturb_bf, num_objects_perturb = segment_BF(im_perturb_bf)
labeled_img_dmso_bf, num_objects_dmso = segment_BF(im_dmso_bf)
dmso_intensities = integrated_intensity(im_dmso_bf, labeled_img_dmso_bf)
perturb_intensities = integrated_intensity(im_perturb_bf, labeled_img_perturb_bf)

# plt.imshow(labeled_img_perturb_bf)
# plt.show()

pickle.dump(dmso_intensities, open('dmso_intensities.pkl', 'wb'))
pickle.dump(perturb_intensities, open('perturb_intensities.pkl', 'wb'))

# # plt.figure()
# # sns.kdeplot(np.array(dmso_intensities),shade=True,color='r')
# # sns.kdeplot(np.array(perturb_intensities),shade=True,color='b')
# # plt.show()
# # ------------------------------------------------
# # ------------------------------------------------


# # ------------------------------------------------
# # Get ND2
# # ------------------------------------------------
# # IMAGE_NAME = 'Plate000_WellA04_Seq0003.nd2'
# # # im_path = IMAGE_DIR + IMAGE_NAME
# # # im_nd2 = get_nd2_dict(im_path)
# # im_nd2 = get_nd2_dict(IMAGE_NAME)


# # # print im_nd2
# # ------------------------------------------------
# # ------------------------------------------------

# # ------------------------------------------------
# # Nuclei Boundaries
# # ------------------------------------------------
# IMAGE_NAME_PERTURB_BF = 'Plate000_WellE09_Seq0056XY1C1.tif'
# im_DAPI_1 = scipy.misc.imread(IMAGE_DIR + '/' + IMAGE_NAME_PERTURB_BF)
# # IMAGE_NAME = 'Plate000_WellA04_Seq0003.nd2'
# # im_nd2 = get_nd2_dict(IMAGE_NAME)
# # im_BF_1 = im_nd2['BF-Cy3'][1]
# im_DAPI_1_labeled, num_objects = segment_DAPI(im_DAPI_1)


# # # print im_BF_1_labeled[0]
# # plt.imshow(im_BF_1_labeled[0])
# # plt.show()

# find_nuclei_boundaries(im_DAPI_1_labeled, im_DAPI_1)
# # ------------------------------------------------
# # ------------------------------------------------




