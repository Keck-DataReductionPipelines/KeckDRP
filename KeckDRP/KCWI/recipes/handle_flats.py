from .. import KcwiConf
from KeckDRP import conf

import numpy as np
import scipy.interpolate as interp
import pylab as pl
import time

def handle_flats(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    if False:
#    if p.in_proctab():
        p.log.warning("Already processed")
    else:
        p.subtract_bias() #CC
        p.subtract_oscan() #CC
        p.trim_oscan() #CC
        p.correct_gain() #CC
        p.remove_badcols() #CC
        p.remove_crs() #CC
        p.rectify_image() #CC
        # update proc table
        p.update_proctab(suffix='int')
        p.write_proctab()
        # output file
        p.write_image(suffix='int')
        p.log.info("flat reduced")
        # how many flats do we have?
        combine_list = p.n_proctab(target_type='FLATLAMP')
        p.log.info("number of flats = %d" % len(combine_list))
        # create master flat
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_FLATS:
            p.image_combine(combine_list, unit=None, suffix='int',
                          in_directory=conf.REDUXDIR, keylog='FLATLIST')
            # output file and update proc table
            p.update_proctab(suffix='flat_stack', newtype='FLAT')
            p.write_image(suffix='flat_stack')
            p.log.info("flat stack produced")
            idl_reference_procedure = p.get_idl_counterpart(target_type='CONTBARS')
            wavemap = p.read_idl_copy(idl_reference_procedure, suffix='wavemap')
            slicemap = p.read_idl_copy(idl_reference_procedure, suffix='slicemap')
            posmap = p.read_idl_copy(idl_reference_procedure, suffix='posmap')
            newflat = frame
            blueslice = 12
            blueleft = 30
            blueright = 40
            p_order = 7
            qblue = np.where((slicemap.data == blueslice) & (posmap.data >= blueleft) & (posmap.data <= blueright))
            xfb = wavemap.data[qblue]
            yfb = newflat.data[qblue]
            s = np.argsort(xfb)
            xfb = xfb[s]
            yfb = yfb[s]
            invar = 1/(1+np.abs(yfb))
            n = 100
            bkpt = np.min(wavemap.data[qblue]) + np.arange(n + 1) * (np.max(wavemap.data[qblue]) -
                                                                     np.min(wavemap.data[qblue])) / n
            bkpty = interp.griddata(xfb, yfb, bkpt, method = 'cubic')
            flat_fit_coeffs = np.polyfit(bkpt, bkpty, p_order)
#            t, c, k = interp.splrep(bkpt, bkpty, k=p_order)
            flat_fit = np.polyval(flat_fit_coeffs, bkpt)
            # plot data and fit
            pl.ion()
            pl.plot(xfb, yfb)
            pl.plot(bkpt, flat_fit)
            pl.xlabel("angstrom")
            pl.ylabel("counts")
            pl.pause(KcwiConf.PLOTPAUSE)
            pl.clf()
            time.sleep(60)
            p.fit_flat() #CC
            p.update_proctab(suffix='master_flat', newtype='MFLAT')
            p.log.info("master flat produced")
        else:
            p.log.info('need %s flats to produce master' % KcwiConf.MINIMUM_NUMBER_OF_FLATS)
        p.write_proctab()
