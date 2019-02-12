from .. import KcwiConf
from KeckDRP import conf


def make_master_dome_flat(p, frame):
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
        p.log.info("dome flat reduced")
        # how many flats do we have?
        combine_list = p.n_proctab(target_type='DOMEFLAT')
        p.log.info("number of flats = %d" % len(combine_list))
        # create master flat
        if len(combine_list) >= 3:
            p.image_combine(combine_list, unit=None, suffix='int',
                            in_directory=conf.REDUXDIR, keylog='DOMELIST')
            # output file and update proc table
            p.update_proctab(suffix='mdimg', newtype='DOME')
            p.write_image(suffix='mdimg')
            p.log.info("dome flat stack produced")
            p.fit_flat()
            p.update_proctab(suffix='mdome', newtype='MDOME')
            p.log.info("master dome flat produced")
        else:
            p.log.info("need 3 dome flats to produce master")
        p.write_proctab()
