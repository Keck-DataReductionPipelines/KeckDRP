from ..core import CcdPrimitives
from ..core import ImgmathPrimitives
from ..core import ProctabPrimitives
from ..core import DevelopmentPrimitives
import os
from .. import conf
from . import KcwiConf

import numpy as np
import scipy as sp
import scipy.interpolate as interp
from scipy.signal import find_peaks
from scipy.signal.windows import boxcar
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interpolate
from scipy.optimize import curve_fit
from scipy.stats import sigmaclip, mode
from skimage import transform as tf
import pickle
from astropy.table import Table
import astropy.io.fits as pf
import matplotlib.pyplot as pl
import time
import math

from astropy.nddata import VarianceUncertainty

import pkg_resources

import KeckDRP

################
# ccdproc usage
# import ccdproc
################


def pascal_shift(coef=None, x0=None):
    """Shift coefficients to a new reference value (X0)

    This should probably go somewhere else, but will be needed here.
    """
    if not coef:
        print("Error, no coefficients for pascal_shift.")
        return None
    if not x0:
        print("Warning, no reference value (x0) supplied")
        return coef
    if len(coef) == 7:
        usecoeff = list(reversed(coef))
        fincoeff = [0.] * 7
    else:
        if len(coef) > 7:
            print("Warning - this routine only handles up to 7 coefficients.")
            usecoeff = list(reversed(coef[0:7]))
            fincoeff = [0.] * len(coef)
        else:
            usecoeff = [0.] * 7
            fincoeff = usecoeff
            for ic, c in enumerate(coef):
                usecoeff[len(coef)-(ic+1)] = coef[ic]
    # get reference values
    x01 = x0
    x02 = x0**2
    x03 = x0**3
    x04 = x0**4
    x05 = x0**5
    x06 = x0**6
    # use Pascal's Triangle to shift coefficients
    fincoeff[0] = usecoeff[0] - usecoeff[1] * x01 + usecoeff[2] * x02 \
        - usecoeff[3] * x03 + usecoeff[4] * x04 - usecoeff[5] * x05 \
        + usecoeff[6] * x06

    fincoeff[1] = usecoeff[1] - 2.0 * usecoeff[2] * x01 \
        + 3.0 * usecoeff[3] * x02 - 4.0 * usecoeff[4] * x03 \
        + 5.0 * usecoeff[5] * x04 - 6.0 * usecoeff[6] * x05

    fincoeff[2] = usecoeff[2] - 3.0 * usecoeff[3] * x01 \
        + 6.0 * usecoeff[4] * x02 - 10.0 * usecoeff[5] * x03 \
        + 15.0 * usecoeff[6] * x04

    fincoeff[3] = usecoeff[3] - 4.0 * usecoeff[4] * x01 \
        + 10.0 * usecoeff[5] * x02 - 20.0 * usecoeff[6] * x03

    fincoeff[4] = usecoeff[4] - 5.0 * usecoeff[5] * x01 \
        + 15.0 * usecoeff[6] * x02

    fincoeff[5] = usecoeff[5] - 6.0 * usecoeff[6] * x01

    fincoeff[6] = usecoeff[6]
    # Trim if needed
    if len(coef) < 7:
        fincoeff = fincoeff[0:len(coef)]
    # Reverse for python
    return list(reversed(fincoeff))
    # END: pascal_shfit()


def gaus(x, a, mu, sigma):
    """Gaussian fitting function"""
    return a * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))


def get_line_window(y, c, thresh=0., verbose=False):
    """Find a window that includes the fwhm of the line"""
    nx = len(y)
    # check edges
    if c < 2 or c > nx - 2:
        if verbose:
            print("input center too close to edge")
        return None, None, 0
    # get initial values
    x0 = c - 2
    x1 = c + 2
    mx = np.nanmax(y[x0:x1+1])
    count = 5
    # check low side
    if x0 - 1 < 0:
        if verbose:
            print("max check: low edge hit")
        return None, None, 0
    while y[x0-1] > mx:
        x0 -= 1
        count += 1
        if x0 - 1 < 0:
            if verbose:
                print("Max check: low edge hit")
            return None, None, 0

    # check high side
    if x1 + 1 > nx:
        if verbose:
            print("max check: high edge hit")
        return None, None, 0
    while y[x1+1] > mx:
        x1 += 1
        count += 1
        if x1 + 1 > nx:
            if verbose:
                print("Max check: high edge hit")
            return None, None, 0
    # adjust starting window to center on max
    cmx = x0 + y[x0:x1+1].argmax()
    x0 = cmx - 2
    x1 = cmx + 2
    mx = np.nanmax(y[x0:x1 + 1])
    # make sure max is high enough
    if mx < thresh:
        return None, None, 0
    #
    # expand until we get to half max
    hmx = mx * 0.5
    #
    # Low index side
    prev = mx
    while y[x0] > hmx:
        if y[x0] > mx or x0 <= 0 or y[x0] > prev:
            if verbose:
                if y[x0] > mx:
                    print("hafmax check: low index err - missed max")
                if x0 <= 0:
                    print("hafmax check: low index err - at edge")
                if y[x0] > prev:
                    print("hafmax check: low index err - wiggly")
            return None, None, 0
        prev = y[x0]
        x0 -= 1
        count += 1
    # High index side
    prev = mx
    while y[x1] > hmx:
        if y[x1] > mx or x1 >= nx or y[x1] > prev:
            if verbose:
                if y[x1] > mx:
                    print("hafmax check: high index err - missed max")
                if x1 >= nx:
                    print("hafmax check: high index err - at edge")
                if y[x1] > prev:
                    print("hafmax check: high index err - wiggly")
            return None, None, 0
        prev = y[x1]
        x1 += 1
        count += 1
    # where did we end up?
    if c < x0 or x1 < c:
        if verbose:
            print("initial position outside final window")
        return None, None, 0

    return x0, x1, count
    # END: get_line_window()


def findpeaks(x, y, wid, sth, ath, pkg=None, verbose=False):
    """Find peaks in spectrum"""
    # derivative
    grad = np.gradient(y)
    # smooth derivative
    win = boxcar(wid)
    d = sp.signal.convolve(grad, win, mode='same') / sum(win)
    # size
    nx = len(x)
    # set up windowing
    if not pkg:
        pkg = wid
    hgrp = int(pkg/2)
    pks = []
    sgs = []
    # loop over spectrum
    # limits to avoid edges given pkg
    for i in np.arange(pkg, (nx - pkg)):
        # find zero crossings
        if np.sign(d[i]) > np.sign(d[i+1]):
            # pass slope threshhold?
            if (d[i] - d[i+1]) > sth * y[i]:
                # pass amplitude threshhold?
                if y[i] > ath or y[i+1] > ath:
                    # get subvectors around peak in window
                    xx = x[(i-hgrp):(i+hgrp+1)]
                    yy = y[(i-hgrp):(i+hgrp+1)]
                    if len(yy) > 3:
                        try:
                            res, _ = curve_fit(gaus, xx, yy,
                                               p0=[100., x[i], 1.])
                            r = abs(x - res[1])
                            t = r.argmin()
                            if abs(i - t) > pkg:
                                if verbose:
                                    print(i, t, x[i], res[1], x[t])
                            else:
                                pks.append(res[1])
                                sgs.append(abs(res[2]))
                        except RuntimeError:
                            continue
    # clean by sigmas
    cpks = []
    sgmd = None
    if len(pks) > 0:
        cln_sgs, low, upp = sigmaclip(sgs, low=3., high=3.)
        for i in range(len(pks)):
            if low < sgs[i] < upp:
                cpks.append(pks[i])
        # sgmn = cln_sgs.mean()
        sgmd = float(np.nanmedian(cln_sgs))
    else:
        print("No peaks found!")
    return cpks, sgmd
    # END: findpeaks()


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives,
                     ProctabPrimitives, DevelopmentPrimitives):

    def __init__(self):
        # KCWI constants
        self.NBARS = 120    # number of bars in continuum bars images
        self.REFBAR = 57    # bar number of reference bar
        self.PIX = 0.0150   # pixel size in mm
        self.FCAM = 305.0   # focal length of camera in mm
        self.GAMMA = 4.0    # mean out-of-plane angle for diffraction (deg)
        self.WAVEFID = 3000.    # Fiducial wavelength for wavelength bins

        # create_unc() variables
        self.readnoise = None       # readnoise (e-)
        # find_bars() variables
        self.refdelx = None         # Reference bar separation in pixels
        # trace_bars() variables
        self.midrow = None          # middle row
        self.midcntr = None         # middle centroids
        self.xi = None
        self.yi = None
        self.xo = None
        self.yo = None
        self.win = None             # sample window
        self.src = None             # source control points
        self.dst = None             # destination control points
        self.barid = None           # control points bar number
        self.slid = None            # control points slice number
        # extract_arcs() variables
        self.arcs = None            # extracted arcs
        # arc_offsets() variables
        self.baroffs = None         # pixel offsets relative to ref bar
        # calc_prelim_disp() variables
        self.prelim_disp = None     # calculated dispersion
        # read_atlas() variables
        self.xvals = None           # pixel values centered on the middle
        self.x0 = None              # middle pixel
        self.reflux = None          # Atlas spectrum
        self.refwave = None         # Altas wavelengths
        self.refdisp = None         # Atlas dispersion
        self.minrow = None          # Lower limit for central fit (px)
        self.maxrow = None          # Upper limit for central fit (px)
        self.offset_wave = None     # atlas-arc offset in Angstroms
        self.offset_pix = None      # atlas-arc offset in pixels
        self.atrespix = None        # atlas matched resolution in atlas px
        # fit_center() variables
        self.centcoeff = []         # Coeffs for central fit of each bar
        # get_atlas_lines() variables
        self.twkcoeff = []          # Coeffs pascal shifted
        self.atminrow = None        # atlas minimum row
        self.atmaxrow = None        # atlas maximum row
        self.atminwave = None       # atlas minimum wavelength (A)
        self.atmaxwave = None       # atlas maximum wavelength (A)
        self.at_wave = None         # atlas wavelength list
        # solve_arcs() variables
        self.fincoeff = []          # Final wavelength solution coeffs
        self.bar_sig = []           # Final sigma of bar wavelength fit
        self.bar_nls = []           # Final number of lines used for fit
        self.xsvals = None          # pixels values starting at zero
        # solve_geom() variables
        self.arc_xpos = None        # arc ref bar line x position list
        self.arc_wave = None        # arc ref var line wavelength list
        self.wavegood0 = None       # Minimum wave at which all slices are good
        self.wavegood1 = None       # Maximum wave at which all slices are good
        self.waveall0 = None        # Minimum wave at which slice data exists
        self.waveall1 = None        # Maximum wave at which slice data exists
        self.wavemid = None         # Middle wavelength
        self.dwout = None           # Output dispersion
        self.wave0out = None        # Output starting wavelength
        self.wave1out = None        # Output ending wavelength
        self.refoutx = None         # Output x positions for bars in cube
        super(KcwiPrimitives, self).__init__()

    @staticmethod
    def kcwi_plot_setup():
        if KcwiConf.INTER >= 1:
            pl.ion()
            fig = pl.figure(num=0, figsize=(17, 6))
            fig.canvas.set_window_title('KCWI DRP')

    def write_image(self, suffix=None):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            if 'BUNIT' in self.frame.header:
                self.frame.unit = self.frame.header['BUNIT']
            outfn = os.path.join(conf.REDUXDIR,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            if not conf.OVERWRITE and os.path.exists(outfn):
                self.log.error("output file exists: %s" % outfn)
            else:
                self.frame.write(outfn, overwrite=conf.OVERWRITE)
                self.log.info("output file: %s" % outfn)

    def write_table(self, table=None, suffix='table', names=None,
                    comment=None, keywords=None):
        if suffix is not None and table is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(conf.REDUXDIR,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            if not conf.OVERWRITE and os.path.exists(outfn):
                self.log.error("output file exists: %s" % outfn)
            else:
                t = Table(table, names=names)
                t.meta['OFNAME'] = origfn
                if comment:
                    t.meta['COMMENT'] = comment
                if keywords:
                    for k, v in keywords.items():
                        t.meta[k] = v
                t.write(outfn, format='fits')
                self.log.info("output file: %s" % outfn)

    def read_table(self, tab=None, indir=None, suffix=None):
        # Set up return table
        retab = None
        # Construct table file name
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
            for f in flist:
                infile = os.path.join(pref, f.split('.')[0] + suff)
                self.log.info("reading table: %s" % infile)
                retab = Table.read(infile, format='fits')
        else:
            self.log.error("No table to read")
        return retab

    def write_geom(self, suffix=None):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(conf.REDUXDIR,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            if not conf.OVERWRITE and os.path.exists(outfn):
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
            self.frame.header['BIASSUB'] = (False,
                                            self.keyword_comments['BIASSUB'])
            self.log.warn('No Master Bias frame found. NO BIAS SUBTRACTION')

    def create_unc(self):
        """Assumes units of image are electron"""
        # start with Poisson noise
        self.frame.uncertainty = VarianceUncertainty(
            self.frame.data, unit='electron^2', copy=True)
        # add readnoise, if known
        if self.readnoise:
            for ia in range(self.frame.namps()):
                sec, rfor = self.parse_imsec(
                    section_key='ATSEC%d' % (ia + 1))
                self.frame.uncertainty.array[
                 sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)] += self.readnoise[ia]
                self.frame.header['BIASRN%d' % (ia + 1)] = self.readnoise[ia]
        else:
            self.log.warn("Readnoise undefined, uncertainty is Poisson only")
        # document variance image creation
        self.frame.header['UNCVAR'] = (True, "has variance image been created?")
        logstr = self.create_unc.__module__ + "." + \
                 self.create_unc.__qualname__
        self.frame.header['HISTORY'] = logstr
        self.log.info(self.create_unc.__qualname__)

    def subtract_dark(self):
        tab = self.n_proctab(target_type='MDARK', nearest=True)
        self.log.info("%d master dark frames found" % len(tab))
        if len(tab) > 0:
            self.img_subtract(tab, suffix='master_dark', indir='redux',
                              keylog='MDFILE')
            self.frame.header['DARKSUB'] = (True,
                                            self.keyword_comments['DARKSUB'])
            logstr = self.subtract_dark.__module__ + "." + \
                     self.subtract_dark.__qualname__
            self.frame.header['HISTORY'] = logstr
            self.log.info(self.subtract_dark.__qualname__)
        else:
            self.frame.header['DARKSUB'] = (False,
                                            self.keyword_comments['DARKSUB'])
            self.log.warn('No Master Dark frame found. NO DARK SUBTRACTION')

    def fit_flat(self):
        tab = self.n_proctab(target_type='FLAT', nearest=True)
        self.log.info("%d flat stacks found" % len(tab))
        # disable this for now
        if len(tab) > 10:
            idl_reference_procedure = self.get_idl_counterpart(
                target_type='CONTBARS')
            wavemap = self.read_idl_copy(idl_reference_procedure,
                                         suffix='wavemap')
            slicemap = self.read_idl_copy(idl_reference_procedure,
                                          suffix='slicemap')
            posmap = self.read_idl_copy(idl_reference_procedure,
                                        suffix='posmap')
            if posmap is None:
                self.log.warn(
                    "No idl reference file found. Stacking is impossible")
                return
            newflat = self.frame
            blueslice = 12
            blueleft = 30
            blueright = 40
            # p_order = 7
            qblue = np.where((slicemap.data == blueslice) &
                             (posmap.data >= blueleft) &
                             (posmap.data <= blueright))
            xfb = wavemap.data[qblue]
            yfb = newflat.data[qblue]
            s = np.argsort(xfb)
            xfb = xfb[s]
            yfb = yfb[s]
            # invar = 1 / (1 + np.abs(yfb))
            n = 100
            bkpt = np.min(wavemap.data[qblue]) + np.arange(n + 1) * \
                (np.max(wavemap.data[qblue]) - np.min(wavemap.data[qblue])) / n
            #          bkpty = interp.griddata(xfb, yfb, bkpt, method = 'cubic')
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
            self.log.info("fit_flat")
    # END: fit_flat()

    def bias_readnoise(self, tab=None, in_directory=None):
        if tab is not None:
            file_list = tab['OFNAME']
            image_numbers = tab['FRAMENO']

            if in_directory is None:
                prefix = '.'
            else:
                prefix = in_directory
            suffix = '.fits'
            # get second and third image in stack
            infil1 = os.path.join(prefix, file_list[1].split('.')[0] + suffix)
            bias1 = KeckDRP.KcwiCCD.read(infil1, unit='adu')
            bias1.data = bias1.data.astype(np.float64)
            infil2 = os.path.join(prefix, file_list[2].split('.')[0] + suffix)
            bias2 = KeckDRP.KcwiCCD.read(infil2, unit='adu')
            bias2.data = bias2.data.astype(np.float64)
            namps = bias1.header['NVIDINP']
            for ia in range(namps):
                # get gain
                gain = bias1.header['GAIN%d' % (ia + 1)]
                # get amp section
                sec, rfor = self.parse_imsec(
                    section_key='DSEC%d' % (ia + 1))
                diff = bias1.data[sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)] - \
                    bias2.data[sec[0]:(sec[1]+1), sec[2]:(sec[3]+1)]
                diff = np.reshape(diff, diff.shape[0]*diff.shape[1]) * \
                    gain / 1.414

                c, upp, low = sigmaclip(diff, low=3.5, high=3.5)
                bias_rn = c.std()
                self.log.info("Amp%d Read noise from bias in e-: %.3f" %
                              ((ia + 1), bias_rn))
                self.frame.header['BIASRN%d' % (ia + 1)] = \
                    (float("%.3f" % bias_rn), "RN in e- from bias")
                if self.frame.inter() >= 1:
                    if self.frame.inter() >= 2:
                        pl.ion()
                    else:
                        pl.ioff()
                    pl.clf()
                    pl.hist(c, bins=50, range=(-12, 12))
                    ylim = pl.gca().get_ylim()
                    pl.plot([c.mean(), c.mean()], ylim, 'g-.')
                    pl.plot([c.mean()-c.std(), c.mean()-c.std()], ylim, 'r-.')
                    pl.plot([c.mean()+c.std(), c.mean()+c.std()], ylim, 'r-.')
                    pl.xlabel("Bias1 - Bias2 (e-/sqrt(2)")
                    pl.ylabel("Density")
                    pl.gca().margins(0)
                    if self.frame.inter() >= 2:
                        input("Next? <cr>: ")
                    else:
                        pl.pause(self.frame.plotpause())
    # END: bias_readnoise()

    def stack_biases(self):

        # get current group id
        if 'GROUPID' in self.frame.header:
            grpid = self.frame.header['GROUPID'].strip()
        else:
            grpid = None
        # how many biases do we have?
        combine_list = self.n_proctab(target_type='BIAS', target_group=grpid)
        self.log.info("number of biases = %d" % len(combine_list))
        # create master bias
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_BIASES:
            self.image_combine(combine_list, keylog='BIASLIST')
            self.bias_readnoise(combine_list)
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
        if 'GROUPID' in self.frame.header:
            grpid = self.frame.header['GROUPID'].strip()
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

        else:
            self.log.info('need %s flats to produce master' %
                          KcwiConf.MINIMUM_NUMBER_OF_FLATS)
        self.write_proctab()

    def subtract_scattered_light(self):
        # keyword of record
        key = 'SCATSUB'
        if self.frame.nasmask():
            self.log.info("NAS Mask in: skipping scattered light subtraction.")
            self.frame.header[key] = (False, self.keyword_comments[key])
        else:

            # Get size of image
            siz = self.frame.data.shape
            # Get x range for scattered light
            x0 = int(siz[1] / 2 - 180 / self.frame.xbinsize())
            x1 = int(siz[1] / 2 + 180 / self.frame.xbinsize())
            # Get y limits
            y0 = 0
            # y1 = int(siz[0] / 2 - 1)
            # y2 = y1 + 1
            y3 = siz[0]
            # print("x limits: %d, %d, y limits: %d, %d" % (x0, x1, y0, y3))
            # Y data values
            yvals = np.nanmedian(self.frame.data[y0:y3, x0:x1], axis=1)
            # X data values
            xvals = np.arange(len(yvals), dtype=np.float)
            # Break points
            nbkpt = int(siz[1]/40.)
            bkpt = xvals[nbkpt:-nbkpt:nbkpt]
            # B-spline fit
            bspl = sp.interpolate.LSQUnivariateSpline(xvals, yvals, bkpt)
            if KcwiConf.INTER >= 1:
                # plot
                pl.ion()
                pl.plot(xvals, yvals, 'ro')
                legend = ["Scat", ]
                xx = np.linspace(0, max(xvals), len(yvals)*5)
                pl.plot(xx, bspl(xx), 'b-')
                legend.append("fit")
                pl.xlabel("y pixel")
                pl.ylabel("e-")
                pl.title("Scat Light img #%d" % (self.frame.header['FRAMENO']))
                pl.legend(legend)
                if KcwiConf.INTER >= 2:
                    input("Next? <cr>: ")
                else:
                    pl.pause(KcwiConf.PLOTPAUSE)
            # Scattered light vector
            scat = bspl(xvals)
            # Subtract scattered light
            self.log.info("Starting scattered light subtraction")
            for ix in range(0, siz[1]):
                self.frame.data[y0:y3, ix] = \
                    self.frame.data[y0:y3, ix] - scat
            self.frame.header[key] = (True, self.keyword_comments[key])

        logstr = self.subtract_scattered_light.__module__ + "." + \
                 self.subtract_scattered_light.__qualname__
        self.frame.header['HISTORY'] = logstr
        self.log.info(self.subtract_scattered_light.__qualname__)
    # END: subtract_scattered_light

    def find_bars(self):
        self.log.info("Finding continuum bars")
        # Do we plot?
        if KcwiConf.INTER >= 1:
            do_plot = True
            pl.ion()
        else:
            do_plot = False
        # initialize
        midcntr = []
        # get image dimensions
        nx = self.frame.data.shape[1]
        ny = self.frame.data.shape[0]
        # get binning
        ybin = self.frame.ybinsize()
        win = int(10 / ybin)
        # select from center rows of image
        midy = int(ny / 2)
        midvec = np.median(self.frame.data[(midy-win):(midy+win+1), :], axis=0)
        # set threshold for peak finding
        midavg = np.average(midvec)
        self.log.info("peak threshold = %f" % midavg)
        # find peaks above threshold
        midpeaks, _ = find_peaks(midvec, height=midavg)
        # do we have the requisite number?
        if len(midpeaks) != self.NBARS:
            self.log.error("Did not find %d peaks: n peaks = %d"
                           % (self.NBARS, len(midpeaks)))
        else:
            self.log.info("found %d bars" % len(midpeaks))
            if do_plot:
                # plot the peak positions
                pl.plot(midvec, '-')
                pl.plot(midpeaks, midvec[midpeaks], 'rx')
                pl.plot([0, nx], [midavg, midavg], '--', color='grey')
                pl.xlabel("CCD X (px)")
                pl.ylabel("e-")
                pl.title("Img %d, Thresh = %.2f" %
                         (self.frame.header['FRAMENO'], midavg))
            # calculate the bar centroids
            for peak in midpeaks:
                xs = list(range(peak-win, peak+win+1))
                ys = midvec[xs] - np.nanmin(midvec[xs])
                xc = np.sum(xs*ys) / np.sum(ys)
                midcntr.append(xc)
                if do_plot:
                    pl.plot([xc, xc], [midavg, midvec[peak]], '-.',
                            color='grey')
            if do_plot:
                pl.plot(midcntr, midvec[midpeaks], 'gx')
                if KcwiConf.INTER >= 2:
                    input("next: ")
                else:
                    pl.pause(self.frame.plotpause())
            self.log.info("Found middle centroids for continuum bars")
        # store peaks
        self.midcntr = midcntr
        # store where we got them
        self.midrow = midy
        self.win = win
        # calculate reference delta x based on refbar
        self.refdelx = 0.
        for ib in range(self.REFBAR-1, self.REFBAR+3):
            self.refdelx += (midcntr[ib] - midcntr[ib-1])
        self.refdelx /= 4.
        self.log.info("Reference delx = %.2f px" % self.refdelx)
        # calculate reference output x values
        # refoutx = midcntr[self.REFBAR] + np.arange(-2, 3) * self.refdelx
        # x0out = int(self.refdelx/2.) + 1
    # END: find_bars()

    def trace_bars(self):
        self.log.info("Tracing continuum bars")
        if KcwiConf.INTER >= 1:
            do_plot = True
            pl.ion()
        else:
            do_plot = False
        if len(self.midcntr) < 1:
            self.log.error("No bars found")
        else:
            # initialize
            samp = int(80 / self.frame.ybinsize())
            win = self.win
            xi = []     # x input
            xo = []     # x output
            yi = []     # y input (and output)
            barid = []  # bar id number
            slid = []   # slice id number
            # loop over bars
            for barn, barx in enumerate(self.midcntr):
                # nearest pixel to bar center
                barxi = int(barx + 0.5)
                # print("bar number %d is at %.3f" % (barn, barx))
                # middle row data
                xi.append(barx)
                xo.append(barx)
                yi.append(self.midrow)
                barid.append(barn)
                slid.append(int(barn/5))
                # trace up
                samy = self.midrow + samp
                done = False
                while samy < (self.frame.data.shape[0] - win) and not done:
                    ys = np.median(
                        self.frame.data[(samy - win):(samy + win + 1),
                                        (barxi - win):(barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn/5))
                    else:
                        done = True
                    samy += samp
                # trace down
                samy = self.midrow - samp
                done = False
                while samy >= win and not done:
                    ys = np.median(
                        self.frame.data[(samy - win):(samy + win + 1),
                                        (barxi - win):(barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 255:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn / 5))
                    else:
                        done = True
                    # disable for now
                    samy -= samp
            # end loop over bars
            # create source and destination coords
            yo = yi
            dst = np.column_stack((xi, yi))
            src = np.column_stack((xo, yo))
            if do_plot:
                # plot them
                pl.clf()
                # pl.ioff()
                pl.plot(xi, yi, 'x', ms=0.5)
                pl.plot(self.midcntr, [self.midrow]*120, 'x', color='red')
                pl.xlabel("CCD X (px)")
                pl.ylabel("CCD Y (px)")
                pl.title("Img %d" % self.frame.header['FRAMENO'])
                if KcwiConf.INTER >= 2:
                    pl.show()
                    input("next: ")
                else:
                    pl.pause(self.frame.plotpause())
            self.write_table(table=[src, dst, barid, slid],
                             names=('src', 'dst', 'barid', 'slid'),
                             suffix='trace',
                             comment=['Source and destination fiducial points',
                                      'Derived from KCWI continuum bars images',
                                      'For defining spatial transformation'],
                             keywords={'MIDROW': (self.midrow,
                                                  "Middle Row of image"),
                                       'WINDOW': (self.win, "Window for bar"),
                                       'REFDELX': (self.refdelx,
                                                   "Reference bar sep in px")})
            if self.frame.saveintims():
                # fit transform
                self.log.info("Fitting spatial control points")
                tform = tf.estimate_transform('polynomial', src, dst, order=3)
                self.log.info("Transforming bars image")
                warped = tf.warp(self.frame.data, tform)
                # write out warped image
                self.frame.data = warped
                self.write_image(suffix='warped')
                self.log.info("Transformed bars produced")
    # END: trace_bars()

    def extract_arcs(self):
        self.log.info("Extracting arc spectra")
        # Find  and read control points from continuum bars
        tab = self.n_proctab(target_type='CONTBARS', nearest=True)
        self.log.info("%d continuum bars frames found" % len(tab))
        trace = self.read_table(tab=tab, indir='redux', suffix='trace')
        self.src = trace['src']  # source control points
        self.dst = trace['dst']  # destination control points
        self.barid = trace['barid']
        self.slid = trace['slid']
        # Get other items
        self.midrow = trace.meta['MIDROW']
        self.win = trace.meta['WINDOW']
        self.refdelx = trace.meta['REFDELX']

        self.log.info("Fitting spatial control points")
        tform = tf.estimate_transform('polynomial', self.src, self.dst, order=3)

        self.log.info("Transforming arc image")
        warped = tf.warp(self.frame.data, tform)
        # Write warped arcs if requested
        if self.frame.saveintims():
            # write out warped image
            self.frame.data = warped
            self.write_image(suffix='warped')
            self.log.info("Transformed arcs produced")
        # extract arcs
        self.log.info("Extracting arcs")
        arcs = []
        for xyi, xy in enumerate(self.src):
            if xy[1] == self.midrow:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped[:, (xi - self.win):(xi + self.win + 1)], axis=1)
                arc = arc - np.nanmin(arc[100:-100])    # avoid ends
                arcs.append(arc)
        # Did we get the correct number of arcs?
        if len(arcs) == self.NBARS:
            self.log.info("Extracted %d arcs" % len(arcs))
            self.arcs = arcs
        else:
            self.log.error("Did not extract %d arcs, extracted %d" %
                           (self.NBARS, len(arcs)))
    # END: extract_arcs()

    def arc_offsets(self):
        self.log.info("Finding inter-bar offsets")
        if self.arcs is not None:
            # Do we plot?
            if KcwiConf.INTER >= 2:
                do_plot = True
                pl.ion()
            else:
                do_plot = False
            # Compare with reference arc
            refarc = self.arcs[self.REFBAR][:]
            # number of cross-correlation samples (avoiding ends)
            nsamp = len(refarc[10:-10])
            # possible offsets
            offar = np.arange(1-nsamp, nsamp)
            # Collect offsets
            offsets = []
            for na, arc in enumerate(self.arcs):
                # Cross-correlate, avoiding junk on the ends
                xcorr = np.correlate(refarc[10:-10], arc[10:-10], mode='full')
                # Calculate offset
                offset = offar[xcorr.argmax()]
                offsets.append(offset)
                self.log.info("Arc %d Slice %d XCorr shift = %d" %
                              (na, int(na/5), offset))
                # display if requested
                if do_plot:
                    pl.clf()
                    pl.plot(refarc, color='green')
                    pl.plot(np.roll(arc, offset), color='red')
                    pl.ylim(bottom=0.)
                    pl.xlabel("CCD y (px)")
                    pl.ylabel("e-")
                    pl.title(self.frame.plotlabel() +
                             " Arc %d Slice %d XCorr, Shift = %d" %
                             (na, int(na/5), offset))
                    pl.show()
                    q = input("<cr> - Next, q to quit: ")
                    if 'Q' in q.upper():
                        do_plot = False
            self.baroffs = offsets
            if self.frame.inter() >= 1:
                if self.frame.inter() >= 2:
                    pl.ion()
                else:
                    pl.ioff()
                pl.clf()
                pl.plot(offsets, 'd')
                ylim = pl.gca().get_ylim()
                for ix in range(1, 24):
                    sx = ix * 5 - 0.5
                    pl.plot([sx, sx], ylim, 'k-.')
                pl.plot([-1, 120], [0., 0.], 'k--')
                pl.xlabel("Bar #")
                pl.ylabel("Offset (px)")
                pl.title(self.frame.plotlabel())
                pl.gca().margins(0)
                if self.frame.inter() >= 2:
                    input("Next? <cr>: ")
                else:
                    pl.pause(self.frame.plotpause())
        else:
            self.log.error("No extracted arcs found")
    # END: arc_offsets()

    def calc_prelim_disp(self):
        # get binning
        ybin = self.frame.ybinsize()
        # 0 - compute alpha
        prelim_alpha = self.frame.grangle() - 13.0 - self.frame.adjang()
        # 1 - compute preliminary angle of diffraction
        prelim_beta = self.frame.camang() - prelim_alpha
        # 2 - compute preliminary dispersion
        prelim_disp = math.cos(prelim_beta/math.degrees(1.)) / \
            self.frame.rho() / self.FCAM * (self.PIX*ybin) * 1.e4
        prelim_disp *= math.cos(self.GAMMA/math.degrees(1.))
        self.log.info("Initial alpha, beta (deg): %.3f, %.3f" %
                      (prelim_alpha, prelim_beta))
        self.log.info("Initial calculated dispersion (A/binned pix): %.3f" %
                      prelim_disp)
        self.prelim_disp = prelim_disp

    def read_atlas(self):
        # What lamp are we using?
        lamp = self.frame.illum()
        atpath = os.path.join(pkg_resources.resource_filename(
            'KeckDRP.KCWI', 'data/'), "%s.fits" % lamp.lower())
        # Does the atlas file exist?
        if os.path.exists(atpath):
            self.log.info("Reading atlas spectrum in: %s" % atpath)
        else:
            self.log.error("Atlas spectrum not found for %s" % atpath)
        # Read the atlas
        ff = pf.open(atpath)
        reflux = ff[0].data
        refdisp = ff[0].header['CDELT1']
        refwav = np.arange(0, len(reflux)) * refdisp + ff[0].header['CRVAL1']
        ff.close()
        # Convolve with appropriate Gaussian
        resolution = self.frame.resolution(refwave=self.frame.cwave())
        atrespix = resolution / refdisp
        self.log.info("Resolution = %.3f Ang, or %.2f Atlas px" % (resolution,
                                                                   atrespix))
        reflux = gaussian_filter1d(reflux, atrespix/2.354)   # convert to sigma
        # Observed arc spectrum
        obsarc = self.arcs[self.REFBAR]
        # Preliminary wavelength solution
        xvals = np.arange(0, len(obsarc)) - int(len(obsarc)/2)
        obswav = xvals * self.prelim_disp + self.frame.cwave()
        # Get central third
        minow = int(len(obsarc)/3)
        maxow = int(2.*len(obsarc)/3)
        # Unless we are low dispersion, then get central 3 5ths
        if 'BL' in self.frame.grating() or 'RL' in self.frame.grating():
            minow = int(len(obsarc)/5)
            maxow = int(4.*len(obsarc)/5)
        minwav = obswav[minow]
        maxwav = obswav[maxow]
        # Get corresponding ref range
        minrw = [i for i, v in enumerate(refwav) if v >= minwav][0]
        maxrw = [i for i, v in enumerate(refwav) if v <= maxwav][-1]
        # Subsample for cross-correlation
        cc_obsarc = obsarc[minow:maxow].copy()
        cc_obswav = obswav[minow:maxow]
        cc_reflux = reflux[minrw:maxrw].copy()
        cc_refwav = refwav[minrw:maxrw]
        # Resample onto reference wavelength scale
        obsint = interpolate.interp1d(cc_obswav, cc_obsarc, kind='cubic',
                                      bounds_error=False,
                                      fill_value='extrapolate'
                                      )
        cc_obsarc = obsint(cc_refwav)
        # Apply cosign bell taper to both
        cc_obsarc *= signal.windows.tukey(len(cc_obsarc),
                                          alpha=self.frame.taperfrac())
        cc_reflux *= signal.windows.tukey(len(cc_reflux),
                                          alpha=self.frame.taperfrac())
        nsamp = len(cc_refwav)
        offar = np.arange(1 - nsamp, nsamp)
        # Cross-correlate
        xcorr = np.correlate(cc_obsarc, cc_reflux, mode='full')
        # Get central region
        x0c = int(len(xcorr)/3)
        x1c = int(2*(len(xcorr)/3))
        xcorr_central = xcorr[x0c:x1c]
        offar_central = offar[x0c:x1c]
        # Calculate offset
        offset_pix = offar_central[xcorr_central.argmax()]
        offset_wav = offset_pix * refdisp
        self.log.info("Initial arc-atlas offset (px, Ang): %d, %.1f" %
                      (offset_pix, offset_wav))
        if self.frame.inter() >= 1:
            if self.frame.inter() >= 2:
                pl.ion()
            else:
                pl.ioff()
            # Plot
            pl.clf()
            pl.plot(offar_central, xcorr_central)
            ylim = pl.gca().get_ylim()
            pl.plot([offset_pix, offset_pix], ylim, 'g-.')
            pl.xlabel("Offset(px)")
            pl.ylabel("X-corr")
            pl.title("Img # %d (%s), Offset = %d px" %
                     (self.frame.header['FRAMENO'], lamp, offset_pix))
            if self.frame.inter() >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.frame.plotpause())
            # Get central wavelength
            cwave = self.frame.cwave()
            # Set up offset tweaking
            q = 'test'
            while q:
                # Plot the two spectra
                pl.clf()
                pl.plot(obswav[minow:maxow] - offset_wav,
                        obsarc[minow:maxow]/np.nanmax(obsarc[minow:maxow]),
                        '-', label="ref bar (%d)" % self.REFBAR)
                pl.plot(refwav[minrw:maxrw],
                        reflux[minrw:maxrw]/np.nanmax(reflux[minrw:maxrw]),
                        'r-', label="Atlas")
                pl.xlim(np.nanmin(obswav[minow:maxow]),
                        np.nanmax(obswav[minow:maxow]))
                ylim = pl.gca().get_ylim()
                pl.plot([cwave, cwave], ylim, 'g-.', label="CWAVE")
                pl.xlabel("Wave(A)")
                pl.ylabel("Rel. Flux")
                pl.title("Img # %d (%s), Offset = %.1f Ang (%d px)" %
                         (self.frame.header['FRAMENO'], lamp,
                          offset_wav, offset_pix))
                pl.legend()
                if self.frame.inter() >= 2:
                    q = input("Enter: <cr> - next, new offset (px): ")
                    if q:
                        try:
                            offset_pix = int(q)
                            offset_wav = offset_pix * refdisp
                        except ValueError:
                            print("Try again")
                else:
                    pl.pause(self.frame.plotpause())
                    q = None
            self.log.info("Final   arc-atlas offset (px, Ang): %d, %.1f" %
                          (offset_pix, offset_wav))
        # Store atlas spectrum
        self.reflux = reflux
        self.refwave = refwav
        self.atrespix = atrespix
        # Store offsets
        self.offset_pix = offset_pix
        self.offset_wave = offset_wav
        # Store reference dispersion
        self.refdisp = refdisp
        # Store central limits
        self.minrow = minow
        self.maxrow = maxow
        # Store x values
        self.xvals = xvals
        self.x0 = int(len(obsarc)/2)
    # END: read_atlas()

    def fit_center(self):
        """ Fit central region

        At this point we have the offsets between bars and the approximate
        offset from the reference bar to the atlas spectrum and the approximate
        dispersion.
        """
        self.log.info("Finding wavelength solution for central region")
        # Are we interactive?
        if KcwiConf.INTER >= 2:
            do_inter = True
            pl.ion()
        else:
            do_inter = False
        # y binning
        ybin = self.frame.ybinsize()
        # let's populate the 0 points vector
        p0 = self.frame.cwave() + np.array(self.baroffs) * self.prelim_disp \
            - self.offset_wave
        # next we are going to brute-force scan around the preliminary
        # dispersion for a better solution. We will wander 5% away from it.
        max_ddisp = 0.05    # fraction
        # we will try nn values
        nn = (int(max_ddisp*abs(self.prelim_disp)/self.refdisp*(
                self.maxrow-self.minrow)/3.0))
        if nn < 10:
            nn = 10
        if nn > 25:
            nn = 25
        self.log.info("N disp. samples: %d" % nn)
        # dispersions to try
        disps = self.prelim_disp * (1.0 + max_ddisp *
                                    (np.arange(0, nn+1) - nn/2.) * 2.0 / nn)
        # containers for bar-specific values
        bardisp = []
        barshift = []
        centwave = []
        centdisp = []

        # values for central fit
        subxvals = self.xvals[self.minrow:self.maxrow]
        # loop over bars
        for b, bs in enumerate(self.arcs):
            # wavelength coefficients
            coeff = [0., 0., 0., 0., 0.]
            # container for maxima, shifts
            maxima = []
            shifts = []
            # get sub spectrum for this bar
            subspec = bs[self.minrow:self.maxrow]
            # now loop over dispersions
            for di, disp in enumerate(disps):
                # populate the coefficients
                coeff[4] = p0[b]
                coeff[3] = disp
                cosbeta = disp / (self.PIX*ybin) * self.frame.rho() * \
                    self.FCAM * 1.e-4
                if cosbeta > 1.:
                    cosbeta = 1.
                beta = math.acos(cosbeta)
                coeff[2] = -(self.PIX * ybin / self.FCAM) ** 2 * \
                    math.sin(beta) / 2. / self.frame.rho() * 1.e4
                coeff[1] = -(self.PIX * ybin / self.FCAM) ** 3 * \
                    math.cos(beta) / 6. / self.frame.rho() * 1.e4
                coeff[0] = (self.PIX * ybin / self.FCAM) ** 4 * \
                    math.sin(beta) / 24. / self.frame.rho() * 1.e4
                # what are the min and max wavelengths to consider?
                wl0 = np.polyval(coeff, self.xvals[self.minrow])
                wl1 = np.polyval(coeff, self.xvals[self.maxrow])
                minwvl = np.nanmin([wl0, wl1])
                maxwvl = np.nanmax([wl0, wl1])
                # where will we need to interpolate to cross-correlate?
                minrw = [i for i, v in enumerate(self.refwave)
                         if v >= minwvl][0]
                maxrw = [i for i, v in enumerate(self.refwave)
                         if v <= maxwvl][-1]
                subrefwvl = self.refwave[minrw:maxrw]
                # need a copy to avoid altering original
                subrefspec = self.reflux[minrw:maxrw].copy()
                # get bell cosine taper to avoid nasty edge effects
                tkwgt = signal.windows.tukey(len(subrefspec),
                                             alpha=self.frame.taperfrac())
                # apply taper to atlas spectrum
                subrefspec *= tkwgt
                # adjust wavelengths
                waves = np.polyval(coeff, subxvals)
                # interpolate the bar spectrum
                obsint = interpolate.interp1d(waves, subspec, kind='cubic',
                                              bounds_error=False,
                                              fill_value='extrapolate')
                intspec = obsint(subrefwvl)
                # apply taper to bar spectrum
                intspec *= tkwgt
                # get a label
                # cross correlate the interpolated spectrum with the atlas spec
                nsamp = len(subrefwvl)
                offar = np.arange(1 - nsamp, nsamp)
                # Cross-correlate
                xcorr = np.correlate(intspec, subrefspec, mode='full')
                # Get central region
                x0c = int(len(xcorr) / 3)
                x1c = int(2 * (len(xcorr) / 3))
                xcorr_central = xcorr[x0c:x1c]
                offar_central = offar[x0c:x1c]
                # Calculate offset
                maxima.append(xcorr_central[xcorr_central.argmax()])
                shifts.append(offar_central[xcorr_central.argmax()])
            # Get interpolations
            int_max = interpolate.interp1d(disps, maxima, kind='cubic',
                                           bounds_error=False,
                                           fill_value='extrapolate')
            int_shift = interpolate.interp1d(disps, shifts, kind='cubic',
                                             bounds_error=False,
                                             fill_value='extrapolate')
            xdisps = np.linspace(min(disps), max(disps), num=nn*100)
            # get peak values
            maxima_res = int_max(xdisps)
            shifts_res = int_shift(xdisps) * self.refdisp
            bardisp.append(xdisps[maxima_res.argmax()])
            barshift.append(shifts_res[maxima_res.argmax()])
            # update coeffs
            coeff[4] = p0[b] - barshift[-1]
            coeff[3] = bardisp[-1]
            cosbeta = coeff[3] / (self.PIX * ybin) * self.frame.rho() * \
                self.FCAM * 1.e-4
            if cosbeta > 1.:
                cosbeta = 1.
            beta = math.acos(cosbeta)
            coeff[2] = -(self.PIX * ybin / self.FCAM) ** 2 * \
                math.sin(beta) / 2. / self.frame.rho() * 1.e4
            coeff[1] = -(self.PIX * ybin / self.FCAM) ** 3 * \
                math.cos(beta) / 6. / self.frame.rho() * 1.e4
            coeff[0] = (self.PIX * ybin / self.FCAM) ** 4 * \
                math.sin(beta) / 24. / self.frame.rho() * 1.e4
            scoeff = pascal_shift(coeff, self.x0)
            self.log.info("Central Fit: Bar#, Cdisp, Coefs: "
                          "%3d  %.4f  %.2f  %.4f  %13.5e %13.5e" %
                          (b, bardisp[-1], scoeff[4], scoeff[3], scoeff[2],
                           scoeff[1]))
            # store central values
            centwave.append(coeff[4])
            centdisp.append(coeff[3])
            # Store results
            self.centcoeff.append(coeff)
            self.twkcoeff.append(scoeff)

            if self.frame.inter() >= 1:
                # plot maxima
                pl.clf()
                pl.plot(disps, maxima, 'r.', label='Data', ms=8)
                pl.plot(xdisps, int_max(xdisps), '-', label='Interp')
                ylim = pl.gca().get_ylim()
                pl.plot([bardisp[-1], bardisp[-1]], ylim, 'g--',
                        label='Peak Disp')
                pl.plot([self.prelim_disp, self.prelim_disp], ylim, 'r-.',
                        label='Calc Disp')
                pl.xlabel("Central Dispersion (Ang/px)")
                pl.ylabel("X-Corr Peak Value")
                pl.title(self.frame.plotlabel() +
                         "Bar %d, Slice %d" % (b, int(b/5)))
                pl.legend()
                pl.gca().margins(0)
                if do_inter:
                    q = input("<cr> - Next, q to quit: ")
                    if 'Q' in q.upper():
                        do_inter = False
                        pl.ioff()
                else:
                    pl.pause(0.01)

        if self.frame.inter() >= 1:
            if self.frame.inter() >= 2:
                pl.ion()
            else:
                pl.ioff()
            # Get central wavelength
            cwave = self.frame.cwave()
            # Plot results
            pl.clf()
            pl.plot(centwave, 'h', label="Data")
            ylim = pl.gca().get_ylim()
            for ix in range(1, 24):
                sx = ix*5 - 0.5
                if ix == 1:
                    pl.plot([sx, sx], ylim, '-.', color='black', label="Slices")
                else:
                    pl.plot([sx, sx], ylim, '-.', color='black')
            pl.xlim([-1, 120])
            pl.plot([-1, 120], [cwave, cwave], 'g-.', label="CWAVE")
            pl.gca().margins(0)
            pl.xlabel("Bar #")
            pl.ylabel("Central Wavelength (A)")
            pl.title(self.frame.plotlabel())
            pl.legend()
            if self.frame.inter() >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.frame.plotpause())
            pl.clf()
            pl.plot(centdisp, 'h', label="Data")
            pl.xlim([-1, 120])
            pl.plot([-1, 120], [self.prelim_disp, self.prelim_disp], 'r-.',
                    label='Calc Disp')
            ylim = pl.gca().get_ylim()
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                if ix == 1:
                    pl.plot([sx, sx], ylim, '-.', color='black', label="Slices")
                else:
                    pl.plot([sx, sx], ylim, '-.', color='black')
            pl.gca().margins(0)
            pl.xlabel("Bar #")
            pl.ylabel("Central Dispersion (A/px)")
            pl.title(self.frame.plotlabel())
            pl.legend()
            if self.frame.inter() >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.frame.plotpause())
    # END: fit_center()

    def get_atlas_lines(self):
        """Get relevant atlas line positions and wavelengths"""
        if KcwiConf.INTER >= 3:
            do_inter = True
            pl.ion()
        else:
            do_inter = False

        # get atlas wavelength range
        #
        # get pixel values (no longer centered in the middle)
        specsz = len(self.arcs[self.REFBAR])
        xvals = np.arange(0, specsz)
        # min, max rows
        minrow = 50
        maxrow = specsz-50
        # wavelength range
        mnwvs = []
        mxwvs = []
        # Get wavelengths
        for b in range(self.NBARS):
            waves = np.polyval(self.twkcoeff[b], xvals)
            mnwvs.append(np.min(waves))
            mxwvs.append(np.max(waves))
        minwav = min(mnwvs) + 10.
        maxwav = max(mxwvs) - 10.
        # Get corresponding atlas range
        minrw = [i for i, v in enumerate(self.refwave) if v >= minwav][0]
        maxrw = [i for i, v in enumerate(self.refwave) if v <= maxwav][-1]
        self.log.info("Min, Max wave (A): %.2f, %.2f" % (minwav, maxwav))
        # store atlas ranges
        self.atminrow = minrw
        self.atmaxrow = maxrw
        self.atminwave = minwav
        self.atmaxwave = maxwav
        # get atlas sub spectrum
        atspec = self.reflux[minrw:maxrw]
        atwave = self.refwave[minrw:maxrw]
        # get reference bar spectrum
        subxvals = xvals[minrow:maxrow]
        subyvals = self.arcs[self.REFBAR][minrow:maxrow].copy()
        subwvals = np.polyval(self.twkcoeff[self.REFBAR], subxvals)
        # smooth subyvals
        win = boxcar(3)
        subyvals = sp.signal.convolve(subyvals, win, mode='same') / sum(win)
        # find good peaks in arc spectrum
        smooth_width = 4                                            # in pixels
        ampl_thresh = 0.
        # slope_thresh = 0.7 * smooth_width / 1000.   # more severe for arc
        # for fitting peaks
        peak_width = int(self.frame.resolution() / self.refdisp)    # in pixels
        if peak_width < 4:
            peak_width = 4
        # slope_thresh = peak_width / 12000.
        slope_thresh = 0.016 / peak_width
        self.log.info("Using a peak_width of %d px, a slope_thresh of %.5f "
                      "a smooth_width of %d and an ampl_thresh of %.3f" %
                      (peak_width, slope_thresh, smooth_width, ampl_thresh))
        init_cent, avwsg = findpeaks(atwave, atspec, smooth_width,
                                     slope_thresh, ampl_thresh, peak_width)
        avwfwhm = avwsg * 2.354
        self.log.info("Found %d peaks with <sig> = %.3f (A), <FWHM> = %.3f (A)"
                      % (len(init_cent), avwsg, avwfwhm))
        if 'BH' in self.frame.grating() or 'BM' in self.frame.grating():
            fwid = avwfwhm
        else:
            fwid = avwfwhm
        # clean near neighbors
        diffs = np.diff(init_cent)
        spec_cent = []
        rej_neigh_w = []
        neigh_fact = 1.25
        for i, w in enumerate(init_cent):
            if i == 0:
                if diffs[i] < avwfwhm * neigh_fact:
                    rej_neigh_w.append(w)
                    continue
            elif i == len(diffs):
                if diffs[i-1] < avwfwhm * neigh_fact:
                    rej_neigh_w.append(w)
                    continue
            else:
                if diffs[i-1] < avwfwhm * neigh_fact or \
                        diffs[i] < avwfwhm * neigh_fact:
                    rej_neigh_w.append(w)
                    continue
            spec_cent.append(w)
        self.log.info("Found %d isolated peaks" % len(spec_cent))
        #
        # generate an atlas line list
        refws = []
        refas = []
        rej_fit_w = []
        rej_par_w = []
        rej_par_a = []
        nrej = 0
        for i, pk in enumerate(spec_cent):

            # Fit Atlas Peak
            line_x = [i for i, v in enumerate(atwave) if v >= pk][0]
            minow, maxow, count = get_line_window(atspec, line_x)
            if count < 5 or not minow or not maxow:
                rej_fit_w.append(pk)
                nrej += 1
                self.log.info("Atlas window rejected for line %.3f" % pk)
                continue
            yvec = atspec[minow:maxow + 1]
            xvec = atwave[minow:maxow + 1]
            try:
                fit, _ = curve_fit(gaus, xvec, yvec, p0=[100., pk, 1.])
            except RuntimeError:
                rej_fit_w.append(pk)
                nrej += 1
                self.log.info("Atlas Gaussian fit rejected for line %.3f" % pk)
                continue
            int_line = interpolate.interp1d(xvec, yvec, kind='cubic',
                                            bounds_error=False,
                                            fill_value='extrapolate')
            x_dense = np.linspace(min(xvec), max(xvec), num=1000)
            # get peak value
            y_dense = int_line(x_dense)
            pki = y_dense.argmax()
            pkw = x_dense[pki]
            # pka = yvec.argmax()
            xoff = abs(pkw - fit[1]) / self.refdisp     # in pixels
            woff = abs(pkw - pk)                        # in Angstroms
            wrat = abs(fit[2]) / fwid                   # can be neg or pos
            if woff > 1. or xoff > 1. or wrat > 1.1:
                rej_par_w.append(pkw)
                rej_par_a.append(y_dense[pki])
                nrej += 1
                self.log.info("Atlas line parameters rejected for line %.3f" %
                              pk)
                self.log.info("woff = %.3f, xoff = %.2f, wrat = %.3f" %
                              (woff, xoff, wrat))
                continue
            refws.append(pkw)
            refas.append(y_dense[pki])
        # store wavelengths
        self.at_wave = refws
        # plot results
        if self.frame.inter() >= 2:
            pl.ion()
        else:
            pl.ioff()
        pl.clf()
        norm_fac = np.nanmax(atspec)
        pl.plot(subwvals, subyvals / np.nanmax(subyvals), label='RefArc')
        xlim = pl.gca().get_xlim()
        ylim = pl.gca().get_ylim()
        pl.plot(atwave, atspec / norm_fac, label='Atlas')
        pl.xlim(xlim)
        # Initial findpeaks list
        plot_first = True
        for iw, w in enumerate(init_cent):
            if plot_first:
                pl.plot([w, w], ylim, 'c--', label='Fpks')
                plot_first = False
            else:
                pl.plot([w, w], ylim, 'c--')
        # Nearby Neighbor rejections
        plot_first = True
        for iw, w in enumerate(rej_neigh_w):
            if plot_first:
                pl.plot([w, w], ylim, 'r--', label='NeighRej')
                plot_first = False
            else:
                pl.plot([w, w], ylim, 'r--')
        # Fit failure rejections
        plot_first = True
        for iw, w in enumerate(rej_fit_w):
            if plot_first:
                pl.plot([w, w], ylim, 'm-.', label='FitRej')
                plot_first = False
            else:
                pl.plot([w, w], ylim, 'm-.')
        # Line Parameter rejections
        plot_first = True
        for iw, w in enumerate(rej_par_w):
            if plot_first:
                pl.plot([w, w], [rej_par_a[iw], rej_par_a[iw]] / norm_fac, 'rd',
                        label='ParRej')
                plot_first = False
            else:
                pl.plot([w, w], [rej_par_a[iw], rej_par_a[iw]] / norm_fac, 'rd')
        # Final Kept list
        plot_first = True
        for iw, w in enumerate(self.at_wave):
            if plot_first:
                pl.plot([w, w], [refas[iw], refas[iw]] / norm_fac, 'kd',
                        label='Kept')
                plot_first = False
            else:
                pl.plot([w, w], [refas[iw], refas[iw]] / norm_fac, 'kd')
        pl.xlabel("Wavelength (A)")
        pl.ylabel("Flux (e-)")
        pl.title(self.frame.plotlabel() + " Ngood = %d, Nrej = %d" %
                 (len(self.at_wave), nrej))
        pl.legend()
        pl.savefig("atlas_lines_%s_%s_%s_%05d.png" %
                   (self.frame.illum(), self.frame.grating(),
                    self.frame.ifuname(), self.frame.header['FRAMENO']))
        self.log.info("Final atlas list has %d lines" % len(self.at_wave))
        if self.frame.inter() >= 2:
            input("Next? <cr>: ")
        else:
            pl.pause(self.frame.plotpause())
    # END: get_atlas_lines()

    def solve_arcs(self):
        """Solve the bar arc wavelengths"""
        if KcwiConf.INTER >= 2:
            master_inter = True
        else:
            master_inter = False
        if KcwiConf.INTER >= 3:
            do_inter = True
            pl.ion()
        else:
            do_inter = False
        verbose = False
        # set thresh for finding lines
        hgt = 50.
        self.log.info("line thresh = %.2f" % hgt)
        # get relevant part of atlas spectrum
        atwave = self.refwave[self.atminrow:self.atmaxrow]
        atspec = self.reflux[self.atminrow:self.atmaxrow]
        # get x values starting at zero pixels
        self.xsvals = np.arange(0, len(self.arcs[self.REFBAR]))
        # loop over arcs and generate a wavelength solution for each
        for ib, b in enumerate(self.arcs):
            # print("")
            # self.log.info("BAR %d" % ib)
            # print("")
            # Starting pascal shifted coeffs from fit_center()
            coeff = self.twkcoeff[ib]
            # get bar wavelengths
            bw = np.polyval(coeff, self.xsvals)
            # smooth spectrum according to slicer
            if 'Small' in self.frame.ifuname():
                bspec = b
            else:
                if 'Large' in self.frame.ifuname():
                    win = boxcar(5)
                else:
                    win = boxcar(3)
                bspec = sp.signal.convolve(b, win, mode='same') / sum(win)
            # spmode = mode(np.round(bspec))
            # spmed = np.nanmedian(bspec)
            # self.log.info("Arc spec median = %.3f, mode = %d" %
            #              (spmed, spmode.mode[0]))
            # store values to fit
            at_wave_dat = []
            arc_pix_dat = []
            rej_wave = []
            nrej = 0
            # loop over lines
            for iw, aw in enumerate(self.at_wave):
                # get window for this line
                try:
                    line_x = [i for i, v in enumerate(bw) if v >= aw][0]
                    minow, maxow, count = get_line_window(bspec, line_x,
                                                          thresh=hgt)
                    if count < 5 or not minow or not maxow:
                        rej_wave.append(aw)
                        nrej += 1
                        if verbose:
                            self.log.info("Arc window rejected for line %.3f"
                                          % aw)
                        continue
                    # check if window no longer contains initial value
                    if minow > line_x > maxow:
                        rej_wave.append(aw)
                        nrej += 1
                        if verbose:
                            self.log.info(
                                "Arc window wandered off for line %.3f" % aw)
                        continue
                    yvec = bspec[minow:maxow+1]
                    xvec = self.xsvals[minow:maxow+1]
                    wvec = bw[minow:maxow+1]
                    max_value = yvec[yvec.argmax()]
                    # Gaussian fit
                    try:
                        fit, _ = curve_fit(gaus, xvec, yvec,
                                           p0=[100., line_x, 1.])
                    except RuntimeError:
                        nrej += 1
                        if verbose:
                            self.log.info(
                                "Arc Gaussian fit rejected for line %.3f" % aw)
                        continue
                    sp_pk_x = fit[1]
                    # Get interpolation
                    int_line = interpolate.interp1d(xvec, yvec, kind='cubic',
                                                    bounds_error=False,
                                                    fill_value='extrapolate')
                    xplot = np.linspace(min(xvec), max(xvec), num=1000)
                    # get peak value
                    plt_line = int_line(xplot)
                    max_index = plt_line.argmax()
                    peak = xplot[max_index]
                    # Calculate centroid
                    cent = np.sum(xvec * yvec) / np.sum(yvec)
                    if abs(cent - peak) > 0.7:
                        rej_wave.append(aw)
                        nrej += 1
                        if verbose:
                            self.log.info("Arc peak - cent offset rejected for "
                                          "line %.3f" % aw)
                        continue
                    # store data
                    arc_pix_dat.append(peak)
                    at_wave_dat.append(aw)
                    # plot, if requested
                    if do_inter:
                        ptitle = "Bar: %d - %3d/%3d: x0, x1, Cent, Wave = " \
                                 "%d, %d, %8.1f, %9.2f" % \
                                 (ib, (iw + 1), len(self.at_wave),
                                  minow, maxow, cent, aw)
                        atx0 = [i for i, v in enumerate(atwave) if v >= min(wvec)][0]
                        atx1 = [i for i, v in enumerate(atwave) if v >= max(wvec)][0]
                        atnorm = np.nanmax(yvec) / np.nanmax(atspec[atx0:atx1])
                        pl.clf()
                        pl.plot(wvec, yvec, 'k--', label='Arc')
                        pl.plot(wvec, yvec, 'r.')
                        ylim = [0, pl.gca().get_ylim()[1]]
                        pl.plot(atwave[atx0:atx1], atspec[atx0:atx1] * atnorm,
                                'g-.', label='Atlas')
                        pl.plot([aw, aw], ylim, 'r-.', label='W in')
                        pl.xlabel("Wavelength (A)")
                        pl.ylabel("Relative Flux")
                        pl.ylim(ylim)
                        pl.title(ptitle)
                        pl.legend()
                        input("next - <cr>: ")
                        pl.clf()
                        pl.plot(xvec, yvec, 'r.', label='Data')
                        pl.plot(xplot, plt_line, label='Interp')
                        ylim = [0, pl.gca().get_ylim()[1]]
                        xlim = pl.gca().get_xlim()
                        pl.plot(xlim, [max_value*0.5, max_value*0.5], 'k--')
                        pl.plot([cent, cent], ylim, 'g--', label='Cntr')
                        pl.plot([line_x, line_x], ylim, 'r-.', label='X in')
                        pl.plot([peak, peak], ylim, 'c-.', label='Peak')
                        pl.plot([sp_pk_x, sp_pk_x], ylim, 'm-.', label='Gpeak')
                        pl.xlabel("CCD Y (px)")
                        pl.ylabel("Flux (DN)")
                        pl.ylim(ylim)
                        pl.title(ptitle)
                        pl.legend()

                        q = input(ptitle + "; <cr> - Next, q to quit: ")
                        if 'Q' in q.upper():
                            do_inter = False
                            pl.ioff()
                except IndexError:
                    if verbose:
                        self.log.info("Atlas line not in observation: %.2f" % aw)
                    rej_wave.append(aw)
                    nrej += 1
                    continue
                except ValueError:
                    if verbose:
                        self.log.info("Interpolation error for line at %.2f" % aw)
                    rej_wave.append(aw)
                    nrej += 1
            self.log.info("Fitting wavelength solution starting with %d "
                          "lines after rejecting %d lines" %
                          (len(arc_pix_dat), nrej))
            # Fit wavelengths
            # Initial fit
            wfit = np.polyfit(arc_pix_dat, at_wave_dat, 4)
            pwfit = np.poly1d(wfit)
            arc_wave_fit = pwfit(arc_pix_dat)
            resid = arc_wave_fit - at_wave_dat
            resid_c, low, upp = sigmaclip(resid, low=3., high=3.)
            wsig = resid_c.std()
            rej_rsd = []
            rej_rsd_wave = []
            # Iteratively remove outliers
            for it in range(4):
                # self.log.info("Iteration %d" % it)
                arc_dat = []
                at_dat = []
                # Trim outliers
                for il, rsd in enumerate(resid):
                    if low < rsd < upp:
                        arc_dat.append(arc_pix_dat[il])
                        at_dat.append(at_wave_dat[il])
                    else:
                        # self.log.info("REJ: %d, %.2f, %.3f" %
                        #              (il, arc_pix_dat[il], at_wave_dat[il]))
                        rej_rsd_wave.append(at_wave_dat[il])
                        rej_rsd.append(rsd)
                # refit
                arc_pix_dat = arc_dat.copy()
                at_wave_dat = at_dat.copy()
                wfit = np.polyfit(arc_pix_dat, at_wave_dat, 4)
                pwfit = np.poly1d(wfit)
                arc_wave_fit = pwfit(arc_pix_dat)
                resid = arc_wave_fit - at_wave_dat
                resid_c, low, upp = sigmaclip(resid, low=3., high=3.)
                wsig = np.nanstd(resid)
            # store results
            # print("")
            self.log.info("Bar %03d, Slice = %02d, RMS = %.3f, N = %d" %
                          (ib, int(ib / 5), wsig, len(arc_pix_dat)))
            # print("")
            self.fincoeff.append(wfit)
            self.bar_sig.append(wsig)
            self.bar_nls.append(len(arc_pix_dat))
            # plot bar fit residuals
            if master_inter:
                pl.ion()
                pl.clf()
                pl.plot(at_wave_dat, resid, 'd', label='Rsd')
                ylim = pl.gca().get_ylim()
                if rej_rsd_wave:
                    pl.plot(rej_rsd_wave, rej_rsd, 'rd', label='Rej')
                pl.xlabel("Wavelength (A)")
                pl.ylabel("Fit - Inp (A)")
                pl.title(self.frame.plotlabel() +
                         " Bar = %03d, Slice = %02d, RMS = %.3f, N = %d" %
                         (ib, int(ib / 5), wsig, len(arc_pix_dat)))
                xlim = [self.atminwave, self.atmaxwave]
                pl.plot(xlim, [0., 0.], '-')
                pl.plot(xlim, [wsig, wsig], '-.', color='gray')
                pl.plot(xlim, [-wsig, -wsig], '-.', color='gray')
                pl.plot([self.frame.cwave(), self.frame.cwave()], ylim, '-.',
                        label='CWAV')
                pl.xlim(xlim)
                pl.ylim(ylim)
                pl.legend()
                input("Next? <cr>: ")

                # overplot atlas and bar using fit wavelengths
                pl.clf()
                bwav = pwfit(self.xsvals)
                pl.plot(bwav, b, label='Arc')
                ylim = pl.gca().get_ylim()
                atnorm = np.nanmax(b) / np.nanmax(atspec)
                pl.plot(atwave, atspec * atnorm, label='Atlas')
                pl.plot([self.frame.cwave(), self.frame.cwave()], ylim, 'm-.',
                        label='CWAV')
                pl.xlim(xlim)
                pl.ylim(ylim)
                pl.xlabel("Wavelength (A)")
                pl.ylabel("Flux")
                pl.title(self.frame.plotlabel() +
                         " Bar = %03d, Slice = %02d, RMS = %.3f, N = %d" %
                         (ib, int(ib/5), wsig, len(arc_pix_dat)))
                leg_first = True
                for w in self.at_wave:
                    if leg_first:
                        pl.plot([w, w], ylim, 'k-.', label='Orig')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'k-.')
                leg_first = True
                for w in arc_wave_fit:
                    if leg_first:
                        pl.plot([w, w], ylim, 'c-.', label='Kept')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'c-.')
                leg_first = True
                for w in rej_rsd_wave:
                    if leg_first:
                        pl.plot([w, w], ylim, 'r-.', label='RejRsd')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'r-.')
                leg_first = True
                for w in rej_wave:
                    if leg_first:
                        pl.plot([w, w], ylim, 'y-.', label='RejFit')
                        leg_first = False
                    else:
                        pl.plot([w, w], ylim, 'y-.')
                pl.legend()
                q = input("Next? <cr>, q - quit: ")
                if 'Q' in q.upper():
                    master_inter = False
                    pl.ioff()
        # Plot final results
        if self.frame.inter() >= 2:
            pl.ion()
        else:
            pl.ioff()
        # Plot fit sigmas
        pl.clf()
        pl.plot(self.bar_sig, 'd', label='RMS')
        xlim = [-1, 120]
        ylim = pl.gca().get_ylim()
        av_bar_sig = float(np.nanmean(self.bar_sig))
        st_bar_sig = float(np.nanstd(self.bar_sig))
        self.log.info("<STD>     = %.3f +- %.3f (A)" % (av_bar_sig, st_bar_sig))
        pl.plot(xlim, [av_bar_sig, av_bar_sig], 'k--')
        pl.plot(xlim, [(av_bar_sig-st_bar_sig), (av_bar_sig-st_bar_sig)], 'k:')
        pl.plot(xlim, [(av_bar_sig+st_bar_sig), (av_bar_sig+st_bar_sig)], 'k:')
        for ix in range(1, 24):
            sx = ix * 5 - 0.5
            pl.plot([sx, sx], ylim, '-.', color='black')
        pl.xlabel("Bar #")
        pl.ylabel("RMS (A)")
        pl.title(self.frame.plotlabel() + " <RMS>: %.3f +- %.3f" % (av_bar_sig,
                                                                    st_bar_sig))
        pl.xlim(xlim)
        pl.gca().margins(0)
        pl.savefig("arc_%05d_resid_%s_%s_%s.png" %
                   (self.frame.header['FRAMENO'], self.frame.illum(),
                    self.frame.grating(), self.frame.ifuname()))
        if self.frame.inter() >= 2:
            input("Next? <cr>: ")
        else:
            pl.pause(self.frame.plotpause())
        # Plot number of lines fit
        pl.clf()
        pl.plot(self.bar_nls, 'd', label='N lines')
        ylim = pl.gca().get_ylim()
        av_bar_nls = float(np.nanmean(self.bar_nls))
        st_bar_nls = float(np.nanstd(self.bar_nls))
        self.log.info("<N Lines> = %.1f +- %.1f" % (av_bar_nls, st_bar_nls))
        pl.plot(xlim, [av_bar_nls, av_bar_nls], 'k--')
        pl.plot(xlim, [(av_bar_nls-st_bar_nls), (av_bar_nls-st_bar_nls)], 'k:')
        pl.plot(xlim, [(av_bar_nls+st_bar_nls), (av_bar_nls+st_bar_nls)], 'k:')
        for ix in range(1, 24):
            sx = ix * 5 - 0.5
            pl.plot([sx, sx], ylim, '-.', color='black')
        pl.xlabel("Bar #")
        pl.ylabel("N Lines")
        pl.title(self.frame.plotlabel() + " <N Lines>: %.3f +- %.3f" %
                 (av_bar_nls, st_bar_nls))
        pl.xlim(xlim)
        pl.gca().margins(0)
        pl.savefig("arc_%05d_nlines_%s_%s_%s.png" %
                   (self.frame.header['FRAMENO'], self.frame.illum(),
                    self.frame.grating(), self.frame.ifuname()))
        if self.frame.inter() >= 2:
            input("Next? <cr>: ")
        else:
            pl.pause(self.frame.plotpause())
        # Plot coefs
        ylabs = ['Ang/px^4', 'Ang/px^3', 'Ang/px^2', 'Ang/px', 'Ang']
        for ic in reversed(range(len(self.fincoeff[0]))):
            pl.clf()
            coef = []
            for c in self.fincoeff:
                coef.append(c[ic])
            pl.plot(coef, 'd')
            ylim = pl.gca().get_ylim()
            for ix in range(1, 24):
                sx = ix * 5 - 0.5
                pl.plot([sx, sx], ylim, '-.', color='black')
            pl.xlabel("Bar #")
            pl.ylabel("Coef %d (%s)" % (ic, ylabs[ic]))
            pl.title(self.frame.plotlabel() + " Coef %d" % ic)
            pl.xlim(xlim)
            pl.gca().margins(0)
            pl.savefig("arc_%05d_coef%d_%s_%s_%s.png" %
                       (self.frame.header['FRAMENO'], ic, self.frame.illum(),
                        self.frame.grating(), self.frame.ifuname()))
            if self.frame.inter() >= 2:
                input("Next? <cr>: ")
            else:
                pl.pause(self.frame.plotpause())
    # END: solve_arcs()

    def solve_geom(self):
        """Solve the overall geometry of the IFU"""
        # Get some geometry constraints
        if self.frame.nasmask():
            goody0 = self.frame.shufrows() + 1
            goody1 = goody0 + self.frame.shufrows()
        else:
            goody0 = 0
            goody1 = max(self.xsvals)
        # Calculate wavelength ranges
        y0wvs = []
        y1wvs = []
        # Get wavelength extremes for each bar
        for fcfs in self.fincoeff:
            y0wvs.append(float(np.polyval(fcfs, goody0)))
            y1wvs.append(float(np.polyval(fcfs, goody1)))
        # Now get ensemble extremes
        y0max = max(y0wvs)
        y0min = min(y0wvs)
        y1max = max(y1wvs)
        y1min = min(y1wvs)
        # Cube trimming wavelengths
        trimw0 = y0min
        trimw1 = y1max
        # Check for negative dispersion
        if trimw0 > trimw1:
            trimw0 = y1min
            trimw1 = y0max
        # Calculate output wavelengths
        self.dwout = self.frame.delta_wave_out()
        ndels = int((trimw0 - self.WAVEFID)/self.dwout)
        self.wave0out = self.WAVEFID + float(ndels) * self.dwout
        ndels = int((trimw1 - self.WAVEFID) / self.dwout)
        self.wave1out = self.WAVEFID + float(ndels) * self.dwout
        self.log.info("WAVE RANGE: %.2f - %.2f" %
                      (self.wave0out, self.wave1out))
        # Calculate wavelength limits
        self.wavegood0 = min([y0max, y1max])
        self.wavegood1 = max([y0min, y1min])
        self.waveall0 = min([y0min, y1min])
        self.waveall1 = max([y0max, y1max])
        self.wavemid = np.average([self.wavegood0, self.wavegood1,
                                   self.waveall0, self.waveall1])
        self.log.info("WAVE  GOOD: %.2f - %.2f" % (self.wavegood0,
                                                   self.wavegood1))
        self.log.info("WAVE   ALL: %.2f - %.2f" % (self.waveall0,
                                                   self.waveall1))
        self.log.info("WAVE   MID: %.2f" % self.wavemid)
        # Start setting up slice transforms
        self.refoutx = np.arange(0, 5) * self.refdelx + int(self.refdelx/2.) + 1
        # Variables for output control points
        srcw = []
        max_srcw = 0
        min_srcw = 4096 / self.frame.ybinsize()
        # Loop over source control points
        for ixy, xy in enumerate(self.src):
            # Calculate y wavelength
            yw = float(np.polyval(self.fincoeff[self.barid[ixy]], xy[1]))
            # Convert to output pixels
            yw = (yw - self.wave0out) / self.dwout
            # Calculate extreme values
            if yw > max_srcw:
                max_srcw = yw
            if yw < min_srcw:
                min_srcw = yw
            srcw.append([xy[0], yw])
        # Use extremes to define output size
        ysize = int(max_srcw + min_srcw + 20 / self.frame.ybinsize())
        xsize = int(max(self.refoutx) + 20 / self.frame.xbinsize())
        self.log.info("Output slices will be %d x %d px" % (xsize, ysize))
        # Now loop over slices and get relevant control points for each slice
        # Output variables
        xl0_out = []
        xl1_out = []
        tform_list = []
        # Loop over 24 slices
        for isl in range(0, 24):
            # Get control points
            xw = []
            yw = []
            xi = []
            yi = []
            # Loop over all control points
            for ixy, xy in enumerate(srcw):
                # Only use the ones for this slice
                if self.slid[ixy] == isl:
                    # Index in to reference output x array
                    ib = self.barid[ixy] % 5
                    # Geometrically corrected control points
                    xw.append(self.refoutx[ib])
                    yw.append(xy[1])
                    # Input control points
                    xi.append(self.dst[ixy][0])
                    yi.append(self.dst[ixy][1])
            # get image limits
            xl0 = int(min(xi) - 24 / self.frame.xbinsize())
            xl1 = int(max(xi) + 16 / self.frame.xbinsize())
            # Store for output
            xl0_out.append(xl0)
            xl1_out.append(xl1)
            self.log.info("Slice %d arc image x limits: %d - %d" %
                          (isl, xl0, xl1))
            # adjust control points
            xit = [x - float(xl0) for x in xi]
            # fit transform
            dst = np.column_stack((xit, yi))
            src = np.column_stack((xw, yw))
            self.log.info("Fitting wavelength and spatial control points")
            tform = tf.estimate_transform('polynomial', src, dst, order=3)
            # Store for output
            tform_list.append(tform)
        # Package geometry data
        geom = {"xsize": xsize, "ysize": ysize,
                "xl0": xl0_out, "xl1": xl1_out,
                "tform": tform_list}
        ofname = self.frame.header['OFNAME']
        outfn = os.path.join(conf.REDUXDIR, ofname.split('.')[0] + '_geom.pkl')
        if os.path.exists(outfn):
            self.log.error("Geometry file already exists: %s" % outfn)
        else:
            with open(outfn, 'wb') as ofile:
                pickle.dump(geom, ofile)
            self.log.info("Geometry written to: %s" % outfn)

    def apply_flat(self):
        self.log.info("apply_flat")

    def subtract_sky(self):
        self.log.info("subtract_sky")

    def make_cube(self):
        self.log.info("Generating data cube")
        # Find and read geometry transformation
        tab = self.n_proctab(target_type='ARCLAMP', nearest=True)
        self.log.info("%d arc frames found" % len(tab))
        ofname = tab['OFNAME'][0]
        geom_file = os.path.join(conf.REDUXDIR,
                                 ofname.split('.')[0] + '_geom.pkl')
        if os.path.exists(geom_file):
            with open(geom_file, 'rb') as ifile:
                geom = pickle.load(ifile)
            xsize = geom['xsize']
            ysize = geom['ysize']
            # out_cube = np.zeros((24, xsize, ysize))
            out_cube = np.zeros((ysize, xsize, 24))
            # set up plots of transformed slices
            pl.clf()
            fig = pl.gcf()
            fig.set_size_inches(5, 12, forward=True)
            # Store original data
            data_img = self.frame.data
            # Loop over 24 slices
            for isl in range(0, 24):
                tform = geom['tform'][isl]
                xl0 = geom['xl0'][isl]
                xl1 = geom['xl1'][isl]
                self.log.info("Transforming image slice %d" % isl)
                slice_img = data_img[:, xl0:xl1]
                pl.clf()
                pl.imshow(slice_img)
                pl.ylim(0, ysize)
                pl.title('raw slice %d' % isl)
                if self.frame.inter() >= 2:
                    input("Next? <cr>: ")
                else:
                    pl.pause(self.frame.plotpause())
                warped = tf.warp(slice_img, tform, order=3,
                                 output_shape=(ysize, xsize))
                print(warped.shape)
                for iy in range(ysize):
                    for ix in range(xsize):
                        out_cube[iy, ix, isl] = warped[iy, ix]
                pl.imshow(warped)
                pl.ylim(0, ysize)
                pl.title('warped slice %d' % isl)
                if self.frame.inter() >= 2:
                    input("Next? <cr>: ")
                else:
                    pl.pause(self.frame.plotpause())
            # write out cube
            self.frame.data = out_cube
            self.write_image(suffix='icube')
            self.log.info("Cube written")
            self.frame.data = data_img
        else:
            self.log.error("Geometry file not found: %s" % geom_file)

    def apply_dar_correction(self):
        self.log.info("apply_dar_correction")

    def flux_calibrate(self):
        self.log.info("flux_calibrate")

    def make_invsensitivity(self):
        self.log.info("make_invsensitivity")
