from .. import KcwiConf


def receive_frame(p, frame):
    # attach frame data
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
    else:
        # update proc table
        p.update_proctab()
        p.write_proctab()
