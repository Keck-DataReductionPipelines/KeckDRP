from ... import conf


def process_darks(p, frame):
    # attach frame data
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab() and conf.OVERWRITE is False:
        p.log.warning("Already processed")
        return

    # process dark frame
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()
    # update proc table
    p.write_image(suffix='int')
    p.update_proctab()
    p.write_proctab()
    p.stack_darks()
