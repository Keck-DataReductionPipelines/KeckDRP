from KeckDRP import PrimitivesBASE
import numpy as np
import scipy as sp
from scipy.signal import find_peaks
import pylab as pl


class SpecPrimitives(PrimitivesBASE):

    def get_bars_pos(self, ypos=None, ncombine=1, threshold=0., pkwin=5):
        # return value
        cntr = []
        # get y limits
        ny = self.frame.data.shape[0]
        # default y position is middle row
        if not ypos:
            ypos = int(ny / 2)
        # check limits
        if ypos - ncombine < 0 or ypos + ncombine + 1 >= ny:
            self.log.error("Y limits exceeded")
        else:
            # extract vector for finding peaks
            vec = np.median(self.frame.data[(ypos - ncombine):
                                            (ypos + (ncombine+1)), :], axis=0)
            peaks, _ = find_peaks(vec, height=threshold)
            # refine positions by taking centroid
            for peak in peaks:
                xs = list(range(peak - pkwin, peak + pkwin + 1))
                ys = vec[xs] - np.nanmin(vec[xs])
                xc = np.sum(xs * ys) / np.sum(ys)
                cntr.append(xc)
