from KCWIPyDRP import conf

def make_master_bias(p, frame):
    # attach frame data
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
    else:
        # update proc table
        p.update_proctab()
        p.write_proctab()
        p.log.info("bias counted")
        # how many biases do we have?
        combine_list = p.n_proctab(targtype='BIAS')
        p.log.info("number of biases = %d" % len(combine_list))
        # create master bias
        if len(combine_list) >= conf.MINIMUM_NUMBER_OF_BIASES:
            p.img_combine(combine_list, keylog='BIASLIST')
            # output file and update proc table
            p.update_proctab(suffix='mbias', newtype='MBIAS')
            p.write_image(suffix='mbias')
            p.log.info("master bias produced")
        else:
            p.log.info(f'need {conf.MINIMUM_NUMBER_OF_BIASES} biases to produce master')
        p.write_proctab()