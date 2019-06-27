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
from skimage import transform as tf
from astropy.table import Table
import pylab as pl
import time
################
# ccdproc usage
# import KeckDRP
# import ccdproc
################


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
                t = Table.read(infile, format='fits')
            return t
        else:
            self.log.error("No table to read")
            return None

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
            ####################################
            # ccdproc usage
            # suff = '_master_bias.fits'
            # pref = 'redux'
            # flist = tab['OFNAME']
            # for f in flist:
            #    infile = os.path.join(pref, f.split('.')[0] + suff)
            #    self.log.info("reading image to subtract: %s" % infile)
            #    master = KeckDRP.KcwiCCD.read(infile, unit='adu')
            # bsub = ccdproc.subtract_bias(self.frame, master)
            # self.set_frame(bsub)
            #####################################
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
            self.img_subtract(tab, suffix='master_dark', indir='redux',
                              keylog='MDFILE')
            self.frame.header['DARKSUB'] = (True,
                                            self.keyword_comments['DARKSUB'])
            logstr = self.subtract_dark.__module__ + "." + \
                     self.subtract_dark.__qualname__
            self.frame.header['HISTORY'] = logstr
            self.log.info(self.subtract_dark.__qualname__)
        else:
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
            self.log.info("fit_flat")

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
            y1 = int(siz[0] / 2 - 1)
            y2 = y1 + 1
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
            # plot
            pl.ion()
            pl.plot(xvals, yvals, 'ro')
            legend = ["Scat",]
            xx = np.linspace(0, max(xvals), len(yvals)*5)
            pl.plot(xx, bspl(xx), 'b-')
            legend.append("fit")
            pl.xlabel("y pixel")
            pl.ylabel("e-")
            pl.title("Scat Light img #%d" % (self.frame.header['FRAMENO']))
            pl.legend(legend)
            if KcwiConf.INTER:
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

    def find_bars(self):
        self.log.info("Finding continuum bars")
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
        if len(midpeaks) != 120:
            self.log.error("Did not find 120 peaks: n peaks = %d"
                           % len(midpeaks))
        else:
            self.log.info("found %d bars" % len(midpeaks))
            # plot the peak positions
            pl.ion()
            pl.plot(midvec, '-')
            pl.plot(midpeaks, midvec[midpeaks], 'rx')
            pl.plot([0, nx], [midavg, midavg], '--', color='grey')
            # calculate the bar centroids
            for peak in midpeaks:
                xs = list(range(peak-win, peak+win+1))
                ys = midvec[xs] - np.nanmin(midvec[xs])
                xc = np.sum(xs*ys) / np.sum(ys)
                pl.plot([xc, xc], [midavg, midvec[peak]], '-.', color='grey')
                midcntr.append(xc)
            pl.plot(midcntr, midvec[midpeaks], 'gx')
            # pl.show()
            pl.pause(self.frame.plotpause())
            self.log.info("Continuum bars centroided")
        self.midcntr = midcntr
        self.midrow = midy
        self.win = win

    def trace_bars(self):
        self.log.info("Tracing continuum bars")
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
                while samy < (self.frame.data.shape[0] - win):
                    ys = np.median(
                        self.frame.data[(samy - win):(samy + win + 1),
                                        (barxi - win):(barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 1500:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn/5))
                    # disable for now
                    if barn == 157:
                        print("bar 57 - xi: %.3f, xo: %.3f, yi: %d" %
                              (xc, barx, samy))
                    samy += samp
                # trace down
                samy = self.midrow - samp
                while samy >= win:
                    ys = np.median(
                        self.frame.data[(samy - win):(samy + win + 1),
                                        (barxi - win):(barxi + win + 1)],
                        axis=0)
                    ys = ys - np.nanmin(ys)
                    xs = list(range(barxi - win, barxi + win + 1))
                    xc = np.sum(xs * ys) / np.sum(ys)
                    if np.nanmax(ys) > 1500:
                        xi.append(xc)
                        xo.append(barx)
                        yi.append(samy)
                        barid.append(barn)
                        slid.append(int(barn / 5))
                    # disable for now
                    if barn == 157:
                        print("bar 57 - xi: %.3f, xo: %.3f, yi: %d" %
                              (xc, barx, samy))
                    samy -= samp
            # end loop over bars
            # create source and destination coords
            yo = yi
            dst = np.column_stack((xi, yi))
            src = np.column_stack((xo, yo))
            # plot them
            pl.clf()
            pl.ion()
            # pl.ioff()
            pl.plot(xi, yi, 'x', ms=0.5)
            pl.plot(self.midcntr, [self.midrow]*120, 'x', color='red')
            # pl.show()
            pl.pause(self.frame.plotpause())
            self.write_table(table=[src, dst], names=('src', 'dst'),
                             suffix='trace',
                             comment=['Source and destination fiducial points',
                                      'Derived from KCWI continuum bars images',
                                      'For defining spatial transformation'],
                             keywords={'MIDROW': self.midrow,
                                       'WINDOW': self.win})
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

    def extract_arcs(self):
        self.log.info("Extracting arc spectra")
        tab = self.n_proctab(target_type='CONTBARS', nearest=True)
        self.log.info("%d continuum bars frames found" % len(tab))
        trace = self.read_table(tab=tab, indir='redux', suffix='trace')
        src = trace['src']
        dst = trace['dst']
        midrow = trace.meta['MIDROW']
        win = trace.meta['WINDOW']
        self.log.info("Fitting spatial control points")
        tform = tf.estimate_transform('polynomial', src, dst, order=3)
        self.log.info("Transforming arc image")
        warped = tf.warp(self.frame.data, tform)
        if self.frame.saveintims():
            # write out warped image
            self.frame.data = warped
            self.write_image(suffix='warped')
            self.log.info("Transformed arcs produced")
        # extract arcs
        self.log.info("Extracting arcs")
        arcs = []
        for xy in src:
            if xy[1] == midrow:
                xi = int(xy[0]+0.5)
                arc = np.median(
                    warped[:, (xi - win):(xi + win + 1)], axis=1)
                arc = arc - np.nanmin(arc)
                pl.clf()
                # pl.ioff()
                pl.ion()
                pl.plot(arc)
                pl.title("Bar pos %.3f" % xy[0])
                # pl.show()
                arcs.append(arc)
                pl.pause(self.frame.plotpause())

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
