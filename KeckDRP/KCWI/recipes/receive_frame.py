from .. import KcwiConf


def receive_frame(p, frame):
    """Process individual frame"""
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return

    # update proc table
    p.update_proctab()
    p.write_proctab()

    p.log.info("frame received")
