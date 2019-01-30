from ..core import CcdPrimitives
from ..core import ImgmathPrimitives
from ..core import ProctabPrimitives
from ..core import DevelopmentPrimitives
import os
from .. import conf


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives,
                     ProctabPrimitives, DevelopmentPrimitives):

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
        tab = self.n_proctab(target_type='MBIAS')
        self.log.info("%d master bias frames found" % len(tab))
        if len(tab) > 0:
            self.img_subtract(tab, suffix='master_bias', indir='redux',
                              keylog='MBFILE')
            self.frame.header['BIASSUB'] = (True,
                                            self.keyword_comments['BIASSUB'])
            logstr = self.subtract_bias.__module__ + "." + \
                     self.subtract_bias.__qualname__
            self.frame.header['HISTORY'] = logstr
            self.log.info(self.subtract_bias.__qualname__)
        else:
            self.log.warn('No Bias frame found. NO BIAS SUBTRACTION')

    def fit_flat(self):
        # idl outline
        # 47 set id of process for ppar comments
        # 50 test: does the given ppar make sense
        # 53 log process and start time
        # 56-59 test: are there flats to fit
        # 62-65 test: is the geometry solved

        # 68-72 read geometry files, test files exist (not needed except to reference next 3)
        # 75-81 read in wavemap, test files exist
        # 84-90 read in slice_image, test files exist
        # 93-99 read in position image, test files exist

        # 102-105 test: directories exist
        # 108-109 list the flat file numbers as integers and count
        # 112 fancy format syntax generator
        # 115 isolate file path (all files must be in same directory it seems)
        # 118 log the number of flats in the stack
        # 121 specify flats are dark subtracted (should be accounted in an earlier primitive in this case)
        # 124-135 test: first file in list exists, or a more primitive version of the first file
        # 138-146 read in first flat, get type and size
        # 149-150 stack images if there is only 1 (just take image)
        # 153-173 stack images if there are only 2
        # 174-218 stack images if there are more than 2
        # 221-231 establish average readnoise
        # 234-235 identify binning mode
        # 238 empty 'ask' string for later
        # 241-248 prep plots if plotting turned on
        # 253-304 corrects a gainfloat issue for a certain configuration. noted as deactivated

        # 309-328 prep variables for actual fitting
            # 328 string = string(list_of_flat_exposure_numbers[0],format = 5 digit integer)
        # 331-449 correct vignetting
            # 333 "get good region for fitting"
            # 340 "get reference slice data"
            # 344 "account for spectral gradient"
            # 349 "fit wavelength slope"
            # 367 "select the points we will fit for the vignetting"
            # 373 "fit vignetted region
            # 382 "plot results"
            # 423 "compute the intersection"
            # 425 "figure out where the correction applies"
            # 427 "apply the correction"
            # 431 "now deal with the intermediate buffer region"
        # 452-461 begin fitting flat
        # 464-532 appended code called "correction for BM where we see a ledge"
        # 534-547 fit the blue end of one slice

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
