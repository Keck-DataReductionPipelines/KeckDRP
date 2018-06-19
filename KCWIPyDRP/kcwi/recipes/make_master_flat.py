from astropy import log


def make_master_flat(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
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
    combine_list = p.n_proctab(targtype='FLATLAMP')
    p.log.info("number of flats = %d" % len(combine_list))
    # create master flat
    if len(combine_list) >= 6:
        p.img_combine(combine_list, suffix='int', indir='redux')
        # output file and update proc table
        p.update_proctab(suffix='mfimg', newtype='FLAT')
        p.write_image(suffix='mfimg')
        p.log.info("master flat produced")
    else:
        p.log.info("need 6 flats to produce master")
    p.write_proctab()
