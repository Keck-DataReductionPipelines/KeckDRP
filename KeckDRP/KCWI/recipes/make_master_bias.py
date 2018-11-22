from .. import KcwiConf

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
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_BIASES:
            master_list = p.n_proctab(targtype='MBIAS')
            if len(master_list) == 0:
                p.img_combine(combine_list, keylog='BIASLIST')
            else:
                p.add_to_master(combine_list, master_list, keylog=None)
            # output file and update proc table
            p.update_proctab(suffix='mbias', newtype='MBIAS')
            p.write_image(suffix='mbias')
            p.log.info("master bias produced")
        else:
            p.log.info('need %s biases to produce master' % (KcwiConf.MINIMUM_NUMBER_OF_BIASES))
        p.write_proctab()