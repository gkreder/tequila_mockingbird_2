import fitness_lib
import pickle
import sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from operator import itemgetter
import matplotlib.colors as colors
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from numpy import ma

# ---------------------------------------------------------------------------------------------------------------------------------------------------
# 
# Plot heatmap with better colorbar
# 
# ---------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------
# Old norm
# --------------------------------------------------------------------
# ## Class to set colermap midpoints to any number
# class MidpointNormalize(colors.Normalize):
#    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#        self.midpoint = midpoint
#        colors.Normalize.__init__(self, vmin, vmax, clip)

#    def __call__(self, value, clip=None):
#        # I'm ignoring masked values and all kinds of edge cases to make a
#        # simple example...
#        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#        return np.ma.masked_array(np.interp(value, x, y))


# --------------------------------------------------------------------
# Iggy Norm
# --------------------------------------------------------------------
class MidpointNormalize(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result


# plt.imshow(p1d1_depth.transpose(), extent=[xmin, xmax, ymin, ymax], cmap='coolwarm', vmin=-0.6, vmax=0.2,
#          norm=MidpointNormalize(midpoint=0.), interpolation='none', aspect='auto')



AMINOTONUMBER_DATA = pickle.load(open('./input_output_files/input/aminotonumber.pkl', 'rb'))
scriptname, replicate = sys.argv
file_name = './input_output_files/output/' + replicate + '_full_data.pkl'


def plot_scores(AA_scores):
    locs = []
    # AAs_keys = AMINOTONUMBER_DATA.keys().
    AA_labels = [tup[0] for tup in sorted(AMINOTONUMBER_DATA.items(), key=itemgetter(1))]
    scores = []

    # numpy.zeros((21,77))
    scores = np.empty((21,75))
    scores[:] = np.NAN
    
    for (loc, AA) in AA_scores:
        # print loc, AA
        if AA != 'WT' and loc > 1 and loc < 77:
            scores[int(AMINOTONUMBER_DATA[AA]), int(loc - 2)] = AA_scores[(loc, AA)]['score']
            locs.append(loc)

    loc_labels = sorted(set(locs))
    # plot_hmap(scores, loc_labels, AA_labels)
    xmin, xmax = 2, 77
    ymin, ymax = 0, 21

    masked_array = np.ma.array(scores, mask=np.isnan(scores))
    # masked_array = scores
    cmap = matplotlib.cm.coolwarm
    cmap.set_bad('black')
    print masked_array


    fig = plt.figure()
    ax = plt.imshow(masked_array, extent=[xmin, xmax, ymin, ymax], cmap=cmap,
            norm=MidpointNormalize(midpoint=0., vmin = -.5, vmax = .25), interpolation='none', aspect='auto')
    # ax = plt.imshow(scores, extent=[xmin, xmax, ymin, ymax], cmap=cmap, vmin=-0.6, vmax=0.2,
        # norm=MidpointNormalize(midpoint=0.), interpolation='none', aspect='auto')

    # Major ticks
    plt.yticks(np.arange(0.5, 21.5, 1), AA_labels[::-1], fontsize=10)
    plt.xticks(np.arange(2.5, 77.5, 1), list(range(xmin, xmax, 1)), fontsize=10)
    # Minor ticks
    plt.axes().set_yticks(np.arange(ymin, ymax, 1), minor=True)
    plt.axes().set_xticks(np.arange(xmin, xmax, 1), minor=True)
    cbar = fig.colorbar(ax)

    plt.axes().spines['bottom'].set_color('w')
    plt.axes().spines['top'].set_color('w') 
    plt.axes().spines['right'].set_color('w')
    plt.axes().spines['left'].set_color('w')
    plt.grid(which='minor', color='w', linestyle='-', linewidth=2)
    plt.grid(which='major', color='w', linestyle='-', linewidth=0)

    plt.tick_params(axis=u'both', which=u'both',length=0)
    
    plt.show()

    # for score in scores:
        # print score

# (AA_scores, total_reads, thrown_out_N_reads, thrown_out_dictionary_reads) = pickle.load(open(file_name, 'rb'))
# plot_scores(AA_scores)

# ---------------------------------------------------------------------------------------------------------------------------------------------------











