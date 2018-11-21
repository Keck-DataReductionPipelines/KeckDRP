from ..core import CcdPrimitives
from ..core import ImgmathPrimitives
from ..core import ProctabPrimitives
import os
from .. import conf



class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives, ProctabPrimitives):

    def __init__(self):
        super(KcwiPrimitives, self).__init__()


    def write_image(self, suffix=None):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(conf.REDUXDIR,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            if not conf.OVERWRITE and os.path.exists(outfn):
                self.log.error("output file exists: %s" % outfn)
            else:
                self.frame.write(outfn, overwrite=conf.OVERWRITE)
                self.log.info("output file: %s" % outfn)

    def write_geom(self, suffix=None):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(conf.REDUXDIR,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            if not self.conf.OVERWRITE and os.path.exists(outfn):
                self.log.error("output file exists: %s" % outfn)
            else:
                # geometry writer goes here
                self.log.info("output file: %s" % outfn)

    def subtract_bias(self):
        tab = self.n_proctab(targtype='MBIAS')
        self.log.info("%d master bias frames found" % len(tab))
        if len(tab)>0:
            self.img_subtract(tab, suffix='mbias', indir='redux', keylog='MBFILE')
            self.frame.header['BIASSUB'] = (True, self.keyword_comments['BIASSUB'])
            logstr = self.subtract_bias.__module__ + "." + \
                 self.subtract_bias.__qualname__
            self.frame.header['HISTORY'] = logstr
            self.log.info(self.subtract_bias.__qualname__)
        else:
            self.log.warn('No Bias frame found. NO BIAS SUBTRACTION')

    def fit_flat(self):
        self.log.info("fit_flat")

    def subtract_scattered_light(self):
        self.log.info("subtract_scattered_light")

    def solve_geom(self):
        self.log.info("solve_geom")

    def apply_flat(self):
        self.log.info("apply_flat")

    def subtract_sky(self):
        self.log.info("subtract_sky")

    def make_cube(self):
        self.log.info("make_cube")

    def apply_dar_correction(self):
        self.log.info("apply_dar_correction")

    def flux_calibrate(self):
        self.log.info("flux_calibrate")

    def make_invsensitivity(self):
        self.log.info("make_invsensitivity")
