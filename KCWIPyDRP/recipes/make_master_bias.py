from astropy import log


def make_master_bias(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    # output file and update proc table
    p.updt_proctab()
    log.info("bias counted")
    # how many biases do we have?
    plist = p.n_proctab(targtype='BIAS')
    log.info("number of biases = %d" % len(plist))
    # create master bias
    if len(plist) >= 7:
        p.img_combine(plist)
        # output file and update proc table
        p.updt_proctab(suffix='mbias', newtype='MBIAS')
        log.info("master bias produced")
    else:
        log.info("need 7 flats to produce master")
    p.write_proctab()
