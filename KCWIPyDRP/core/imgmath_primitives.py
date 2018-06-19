import os
from ..kcwi.kcwi_objects import KcwiCCD
import ccdproc

class ImgmathPrimitives:

    def __init__(self):
        super(ImgmathPrimitives, self).__init__()

    def img_combine(self, tab=None, ctype='bias', indir=None, suffix=None,
                    method='average', unit='adu'):
        if tab is not None:
            flist = tab['OFNAME']

            if indir is None:
                pref = '.'
            else:
                pref = indir

            if suffix is None:
                suff = '.fits'
            else:
                suff = '_' + suffix + '.fits'

            stack = []
            for f in flist:
                infile = os.path.join(pref, f.split('.')[0] + suff)
                self.log.info("reading image: %s" % infile)
                stack.append(KcwiCCD.read(infile, unit=unit))
            if 'bias' in ctype:
                self.set_frame(ccdproc.combine(stack, method=method,
                                               sigma_clip=True,
                                               sigma_clip_low_thresh=None,
                                               sigma_clip_high_thresh=2.0))
            else:
                self.set_frame(ccdproc.combine(stack, method=method))
            self.log.info("img_combine %s using %s" % (type, method))
        else:
            self.log.info("something went wrong with img_combine")

    def img_subtract(self, tab=None, indir=None, suffix=None, unit='adu'):
        if tab is not None:
            flist = tab['OFNAME']

            if indir is None:
                pref = '.'
            else:
                pref = indir

            if suffix is None:
                suff = '.fits'
            else:
                suff = '_' + suffix + '.fits'

            subtrahend = None
            for f in flist:
                infile = os.path.join(pref, f.split('.')[0] + suff)
                self.log.info("reading image to subtract: %s" % infile)
                subtrahend = KcwiCCD.read(infile, unit=unit)

            if subtrahend is not None:
                result = self.frame.subtract(subtrahend)
                result.meta = self.frame.meta

            self.set_frame(result)
            self.log.info("img_subtract")
        else:
            self.log.info("something went wrong with img_subtract")

    def img_divide(self):
        self.log.info("img_divide")
