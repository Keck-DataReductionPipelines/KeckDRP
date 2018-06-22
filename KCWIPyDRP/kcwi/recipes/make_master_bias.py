from astropy import log


def make_master_bias(p, frame):
    # attach frame data
    p.set_frame(frame)
    p.read_proctab()
    # update proc table
    p.update_proctab()
    p.write_proctab()
    p.log.info("bias counted")
    # how many biases do we have?
    combine_list = p.n_proctab(targtype='BIAS')
    p.log.info("number of biases = %d" % len(combine_list))
    # create master bias
    if len(combine_list) >= 7:
        p.img_combine(combine_list, keylog='BIASLIST')
        # output file and update proc table
        p.update_proctab(suffix='mbias', newtype='MBIAS')
        p.write_image(suffix='mbias')
        p.log.info("master bias produced")
    else:
        p.log.info("need 7 biases to produce master")
    p.write_proctab()
