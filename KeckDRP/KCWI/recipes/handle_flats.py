from .. import KcwiConf
from KeckDRP import conf

def handle_flats(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
    else:
        p.subtract_bias()
        p.subtract_oscan()
        p.trim_oscan()
        p.correct_gain()
        p.remove_badcols()
        p.remove_crs()
        p.rectify_image()
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
            p.update_proctab(suffix='mfimg', newtype='FLAT')
            p.write_image(suffix='mfimg')
            p.log.info("flat stack produced")
            p.fit_flat()
            p.update_proctab(suffix='mflat', newtype='MFLAT')
            p.log.info("master flat produced")
        else:
            p.log.info('need %s flats to produce master' % KcwiConf.MINIMUM_NUMBER_OF_FLATS)
        p.write_proctab()
