from .. import KcwiConf


def process_contbars(p, frame):
    """process a continuum bars image and output the bar trace info"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    # reduce contbars
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.rectify_image()

    # write image
    p.write_image(suffix='int')
    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

    p.log.info("contbars reduced")

    # trace contbars
    p.find_bars()
    p.trace_bars()

    p.log.info("spatial geometry traced")
