from astropy import log


def make_master_bias(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    # output file and update proc table
    p.update_proctab()
    p.generate_output_image()
    log.info("bias counted")
    # how many biases do we have?
    combine_list = p.n_proctab(targtype='BIAS')
    log.info("number of biases = %d" % len(combine_list))
    # create master bias
    if len(combine_list) >= 7:
        p.img_combine(combine_list)
        # output file and update proc table
        p.update_proctab(suffix='mbias', newtype='MBIAS')
        p.generate_output_image(suffix='mbias')
        log.info("master bias produced")
    else:
        log.info("need 7 flats to produce master")
    p.write_proctab()
