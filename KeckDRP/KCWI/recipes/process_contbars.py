from .. import KcwiConf


def process_contbars(p, frame):
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.rectify_image()

    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

    # write image
    p.write_image(suffix='int')
    p.log.info("contbars reduced")

    p.find_bars()
    p.trace_bars()

