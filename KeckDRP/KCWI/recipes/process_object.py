

def process_object(p, frame):
    """Process a science image and output a calibrated data cube"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    # basic CCD reduction
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.create_unc()
    p.rectify_image()

    # write image
    p.write_image(suffix='int')
    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

    p.log.info("science frame reduced")

    # Dark/scattered light
    p.subtract_dark()
    p.subtract_scattered_light()

    # write image
    p.write_image(suffix='intd')
    # update proc table
    p.update_proctab(suffix='intd')
    p.write_proctab()

    # Flat field reduction
    p.apply_flat()

    # write image
    # p.write_image(suffix='intf')
    # update proc table
    # p.update_proctab(suffix='intf')
    # p.write_proctab()

    # Sky subtraction
    p.subtract_sky()

    # write image
    # p.write_image(suffix='intk')
    # update proc table
    # p.update_proctab(suffix='intk')
    # p.write_proctab()

    # 3-D reductions
    p.make_cube()

    # write image
    p.write_image(suffix='icube')
    # update proc table
    p.update_proctab(suffix='icube')
    p.write_proctab()

    # DAR correction
    # p.apply_dar_correction()

    # write image
    # p.write_image(suffix='icubed')
    # update proc table
    # p.update_proctab(suffix='icubed')
    # p.write_proctab()

    # Flux calibration
    # p.flux_calibrate()

    # write image
    # p.write_image(suffix='icubes')
    # update proc table
    # p.update_proctab(suffix='icubes')
    # p.write_proctab()

    p.log.info("calibrated science cube generated")
