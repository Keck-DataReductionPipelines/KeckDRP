from ... import conf
from KeckDRP.core.decorators import repeat_processing

@repeat_processing
def process_biases(p, frame):
    # attach frame data
    p.set_frame(frame)
    #p.read_proctab()
    #if p.in_proctab() and conf.OVERWRITE is False:
    #    p.log.warning("Already processed")
    #    return

    # update proc table
    p.update_proctab()
    p.write_proctab()
    p.stack_biases()
