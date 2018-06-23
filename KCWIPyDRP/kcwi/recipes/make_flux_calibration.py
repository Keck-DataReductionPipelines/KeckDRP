

def make_flux_calibration(p, frame):
    # do reduction for a standard star frame
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
    else:
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
        p.make_invsensitivity()
        p.flux_calibrate()
        # this should eventually be suffix='icube'
        p.write_image(suffix='int')
        # update proc table
        p.update_proctab(suffix='int')
        p.write_proctab()
        p.log.info("science frame reduced")

