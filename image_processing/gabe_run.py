from gabe import *
import scipy

IMAGE_NAME_PERTURB_BF = 'Plate000_WellE09_Seq0056XY1C1.tif'
im_DAPI_1 = scipy.misc.imread(IMAGE_DIR + '/' + IMAGE_NAME_PERTURB_BF)
im_DAPI_1_labeled, num_objects = segment_DAPI(im_DAPI_1)

plt.imshow(im_DAPI_1_labeled)
plt.show()