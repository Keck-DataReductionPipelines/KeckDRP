

def process_standard(p, frame):
    """Process standard star observation"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    # reduce standard
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.rectify_image()
    p.apply_flat()
    p.subtract_sky()
    p.make_cube()
    p.apply_dar_correction()

    # write out reduced image
    p.write_image(suffix='int')

    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

    p.log.info("standard star reduced")

    # create calibration
    p.make_invsensitivity()
    p.flux_calibrate()

