from .. import KcwiConf

def handle_biases(p, frame):
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
        combine_list = p.n_proctab(target_type='BIAS')
        p.log.info("number of biases = %d" % len(combine_list))
        # create master bias
        if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_BIASES:
            p.image_combine(combine_list, keylog='BIASLIST')
            # output file and update proc table
            p.update_proctab(suffix='master_bias', newtype='MBIAS')
            p.write_image(suffix='master_bias')
            p.log.info("master bias produced")
        else:
            p.log.info('need %s biases to produce master' % KcwiConf.MINIMUM_NUMBER_OF_BIASES)
        p.write_proctab()
