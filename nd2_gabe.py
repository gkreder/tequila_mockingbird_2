#!/usr/bin/env python2
import nd2reader
from sys import argv
import pickle
import numpy as np

script, nd2 = argv

# nd2 = nd2reader.Nd2(str(nd2))

# dapi = []
# fitc = []
# fitclong = []
# cy3 = []
# bfcy3 = []

# for image in nd2.select(channels = 'DAPI'):
#     dapi.append(image)

# for image in nd2.select(channels = 'FITC'):
#     fitc.append(image)

# for image in nd2.select(channels = 'FITClong'):
#     fitclong.append(image)

# for image in nd2.select(channels = 'CY3'):
#     cy3.append(image)

# for image in nd2.select(channels = 'BF-Cy3'):
#     bfcy3.append(image)

im_name_pickle = str(nd2).replace('.nd2', '.pkl')
nd2 = nd2reader.Nd2(str(nd2))

img_dict = {'DAPI' : {},
			'FITC' : {},
			'FITClong' : {},
			'CY3' : {},
			'BF-Cy3' : {}}

# nd2 = nd2reader.Nd2(IMAGE_DIR+'/'+img_path)


# nd2 = nd2reader.Nd2(img_path)

# print nd2.fields_of_view

for i in range(len(nd2)):
	channel = nd2[i].channel
	fov = nd2[i].field_of_view
	# img_dict[channel][fov] = nd2[i]
	img_dict[channel][fov] = nd2[i].astype(np.uint16) 
	# / np.amax(nd2[i]).astype(np.uint16)

pickle.dump(img_dict, open(im_name_pickle, 'w'))

