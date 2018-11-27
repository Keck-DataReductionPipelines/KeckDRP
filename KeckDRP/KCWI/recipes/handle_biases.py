from .. import KcwiConf

def handle_biases(p, frame):  #CC needs to know frame is a bias in advance!
    # attach frame data
    p.set_frame(frame)
    p.log.info("bias counted")
    # how many biases do we have?
    combine_list = p.n_proctab(targtype='BIAS')
    p.log.info("number of biases = %d" % len(combine_list))
    # create master bias
    if len(combine_list) >= KcwiConf.MINIMUM_NUMBER_OF_BIASES:
        p.img_combine(combine_list, keylog='BIASLIST') #CC if this primitive is unique to biases, let's call it bias_combine.  if it comes from imgmath, let's make a bias_combine in the primitive file using imgcombine.  this will make a good template for future combines and simplify the call.
        # output file and update proc table
        p.update_proctab(suffix='master_bias', newtype='MBIAS')
        p.write_image(suffix='master_bias')
        p.log.info("master bias produced")
    else:
        p.log.info(f'need {KcwiConf.MINIMUM_NUMBER_OF_BIASES} biases to produce master')
    p.write_proctab()
