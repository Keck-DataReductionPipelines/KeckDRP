

def make_master_geom(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
    else:
        p.trim_oscan()
        p.correct_gain()
        p.remove_badcols()
        p.rectify_image()
        # update proc table
        p.update_proctab(suffix='int')
        p.write_proctab()
        # output file
        p.write_image(suffix='int')
        p.log.info("flat reduced")
        # how many geom files do we have?
        arc_list = p.n_proctab(targtype='ARCLAMP')
        p.log.info("number of arcs = %d" % len(arc_list))
        cbars_list = p.n_proctab(targtype='CONTBARS')
        p.log.info("number of cbars = %d" % len(cbars_list))
        # create master geom
        if len(arc_list) >= 2 and len(cbars_list) >= 1:
            p.solve_geom()
            # output file and update proc table
            p.update_proctab(suffix='geom', newtype='MGEOM')
            p.write_geom(suffix='geom')
            p.log.info("master geom produced")
        else:
            p.log.info("need 2 arcs and 1 cbars to produce master")
        p.write_proctab()
