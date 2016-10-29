from sys import argv

script, pkl_filename = argv

import numpy as numpy
import pickle as pkl


differences_array = numpy.load(pkl_filename, 'rb')
#print differences_array
#for i in differences_array:
	#i >= 0.35
	#positive_fitness.append(i)

	#pri
 {k: v for k, v in differences_array.iteritems() if v > 0.35}
 positive fitness =
print positive_fitness