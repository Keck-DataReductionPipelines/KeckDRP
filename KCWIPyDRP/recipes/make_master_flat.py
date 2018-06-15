from astropy import log


def make_master_flat(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()
    # output file and update proc table
    p.updt_proctab(suffix='int')
    log.info("flat reduced")
    # how many flats do we have?
    plist = p.n_proctab(targtype='FLATLAMP')
    log.info("number of flats = %d" % len(plist))
    # create master flat
    if len(plist) >= 6:
        p.img_combine(plist, suffix='int')
        # output file and update proc table
        p.updt_proctab(suffix='mfimg', newtype='FLAT')
        log.info("master flat produced")
    else:
        log.info("need 6 flats to produce master")
    p.write_proctab()
