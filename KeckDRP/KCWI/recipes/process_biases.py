from ... import conf
from KeckDRP.core.decorators import repeat_processing


@repeat_processing
def process_biases(p, frame):
    """Process bias frame, creating master when enough have been taken"""
    p.kcwi_plot_setup()
    p.set_frame(frame)

    # update proc table
    p.update_proctab(suffix='RAW')
    p.write_proctab()

    p.log.info("bias processed")

    # make master bias
    p.stack_biases()
