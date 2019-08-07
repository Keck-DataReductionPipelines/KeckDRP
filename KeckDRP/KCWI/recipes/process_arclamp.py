from .. import KcwiConf


def process_arclamp(p, frame):
    """Process arc lamp images and output geometry solution"""
    p.kcwi_plot_setup()
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    # reduce arc
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
    p.log.info("arclamp reduced")

    # solve geometry
    p.extract_arcs()
    p.arc_offsets()
    p.calc_prelim_disp()
    p.read_atlas()
    p.fit_center()
    p.get_atlas_lines()
    p.solve_geom()
