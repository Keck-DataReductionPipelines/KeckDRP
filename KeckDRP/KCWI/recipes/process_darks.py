from ... import conf


def process_darks(p, frame):
    """Process dark frame, creating master when enough have been taken"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return

    # reduce dark frame
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()
    # write out reduced image
    p.write_image(suffix='int')

    # update proc table
    p.update_proctab()
    p.write_proctab()

    p.log.info("dark reduced")

    # create master dark
    p.stack_darks()
