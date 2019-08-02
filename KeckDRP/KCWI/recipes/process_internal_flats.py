from ... import conf


def process_internal_flats(p, frame):
    """Process internal flat, create master when enough have been taken"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab() and conf.OVERWRITE is False:
        p.log.warning("Already processed")
        return
    # reduce internal flat
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()
    p.subtract_dark()
    p.subtract_scattered_light()

    # output file
    p.write_image(suffix='int')

    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

    p.log.info("flat reduced")

    # create master flat
    p.stack_internal_flats()
    p.fit_flat()


