from .. import KcwiConf
from KeckDRP import conf


def process_dome(p, frame):
    """Process dome flat creating master when enough have been taken"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    # reduce dome flat
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()

    # output file
    p.write_image(suffix='int')

    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

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
