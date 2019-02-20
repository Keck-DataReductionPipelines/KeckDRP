from ..core import CcdPrimitives
from ..core import ImgmathPrimitives
from ..core import ProctabPrimitives
from ..core import DevelopmentPrimitives
import os
from .. import conf
from . import KcwiConf

import numpy as np
import scipy.interpolate as interp
import pylab as pl
import time


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
        tab = self.n_proctab(target_type='MBIAS', nearest=True)
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
            self.log.warn('No Master Bias frame found. NO BIAS SUBTRACTION')

    def subtract_dark(self):
        tab = self.n_proctab(target_type='MDARK', nearest=True)
        self.log.info("%d master dark frames found" % len(tab))
        if len(tab) > 0:
            self.img_subtract(tab, suffix='master_bias', indir='redux',
                              keylog='MDFILE')
            self.frame.header['DARKSUB'] = (True,
                                            self.keyword_comments['DARKSUB'])
            logstr = self.subtract_bias.__module__ + "." + \
                     self.subtract_bias.__qualname__
            self.frame.header['HISTORY'] = logstr
            self.log.info(self.subtract_bias.__qualname__)
        else:
            self.log.warn('No Master Dark frame found. NO DARK SUBTRACTION')

    def fit_flat(self):
        pass

        self.log.info("fit_flat")

    def stack_biases(self):
        # get current group id
        if 'GRPID' in self.frame.header:
            grpid = self.frame.header['GRPID'].strip()
        else:
            grpid = None
        # how many biases do we have?
        combine_list = self.n_proctab(target_type='BIAS', target_group=grpid)
        self.log.info("number of biases = %d" % len(combine_list))
        # create master bias
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_BIASES:
            self.image_combine(combine_list, keylog='BIASLIST')
            # output file and update proc table
            self.update_proctab(suffix='master_bias', newtype='MBIAS')
            self.write_image(suffix='master_bias')
            self.log.info("master bias produced")
        else:
            self.log.info('need %s biases to produce master' %
                          KcwiConf.MINIMUM_NUMBER_OF_BIASES)
        self.write_proctab()

    def stack_darks(self):
        # get current group id
        if 'GRPID' in self.frame.header:
            grpid = self.frame.header['GRPID'].strip()
        else:
            grpid = None
        # how many darks do we have?
        combine_list = self.n_proctab(target_type='DARK', target_group=grpid)
        self.log.info("number of darks = %d" % len(combine_list))
        # create master bias
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_DARKS:
            self.image_combine(combine_list, keylog='DARKLIST',
                               in_directory='redux', suffix='int')
            # output file and update proc table
            self.update_proctab(suffix='master_dark', newtype='MDARK')
            self.write_image(suffix='master_dark')
            self.log.info("master dark produced")
        else:
            self.log.info('need %s darks to produce master' %
                          KcwiConf.MINIMUM_NUMBER_OF_DARKS)
        self.write_proctab()

    def stack_internal_flats(self):
        # how many flats do we have?
        combine_list = self.n_proctab(target_type='FLATLAMP')
        self.log.info("number of flats = %d" % len(combine_list))
        # create master flat
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_FLATS:
            self.image_combine(combine_list, unit=None, suffix='int',
                               in_directory=conf.REDUXDIR, keylog='FLATLIST')
            # output file and update proc table
            self.update_proctab(suffix='flat_stack', newtype='FLAT')
            self.write_image(suffix='flat_stack')
            self.log.info("flat stack produced")
            idl_reference_procedure = self.get_idl_counterpart(
                target_type='CONTBARS')
            wavemap = self.read_idl_copy(idl_reference_procedure,
                                         suffix='wavemap')
            slicemap = self.read_idl_copy(idl_reference_procedure,
                                          suffix='slicemap')
            posmap = self.read_idl_copy(idl_reference_procedure,
                                        suffix='posmap')
            if posmap is None:
                self.log.warn("No idl reference file found. Stacking is impossible")
                return
            newflat = self.frame
            blueslice = 12
            blueleft = 30
            blueright = 40
            p_order = 7
            qblue = np.where((slicemap.data == blueslice) &
                             (posmap.data >= blueleft) &
                             (posmap.data <= blueright))
            xfb = wavemap.data[qblue]
            yfb = newflat.data[qblue]
            s = np.argsort(xfb)
            xfb = xfb[s]
            yfb = yfb[s]
            invar = 1 / (1 + np.abs(yfb))
            n = 100
            bkpt = np.min(wavemap.data[qblue]) + np.arange(n + 1) * \
                (np.max(wavemap.data[qblue]) - np.min(wavemap.data[qblue])) / n
            #            bkpty = interp.griddata(xfb, yfb, bkpt, method = 'cubic')
            bkpty = interp.griddata(xfb, yfb, bkpt)
            t, c, k = interp.splrep(bkpt, bkpty, k=3)
            #            flat_fit_coeffs = np.polyfit(bkpt, bkpty, p_order)
            #            t, c, k = interp.splrep(bkpt, bkpty, k=p_order)
            #            flat_fit = np.polyval(flat_fit_coeffs, bkpt)
            spline = interp.BSpline(t, c, k, extrapolate=False)
            # plot data and fit
            pl.ion()
            pl.plot(xfb, yfb)
            pl.plot(bkpt, spline(bkpt))
            pl.xlabel("angstrom")
            pl.ylabel("counts")
            pl.pause(KcwiConf.PLOTPAUSE)
            pl.clf()
            time.sleep(15)
            self.fit_flat()  # CC
            self.update_proctab(suffix='master_flat', newtype='MFLAT')
            self.log.info("master flat produced")
        else:
            self.log.info('need %s flats to produce master' %
                          KcwiConf.MINIMUM_NUMBER_OF_FLATS)
        self.write_proctab()

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
