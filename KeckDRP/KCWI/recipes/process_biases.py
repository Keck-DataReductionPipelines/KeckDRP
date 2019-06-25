from ... import conf
from KeckDRP.core.decorators import repeat_processing


@repeat_processing
def process_biases(p, frame):
    p.set_frame(frame)
    p.update_proctab(suffix='RAW')
    p.write_proctab()
    p.stack_biases()
