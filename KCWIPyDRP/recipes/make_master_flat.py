from astropy import log
import ccdproc


def make_master_flat(p, frame):
    p.set_frame(frame)
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()
    # p.put_proc()
    # nflat, flist = p.get_proc(type='FLAT')
    #if nflat > 5:
    #    c = p.img_combine(flist)
    #    m = p.fit_flat(c)
    #    p.output_master(m, master_type="FLAT")
    #    p.put_proc(m)
    #    log.info("make_master_flat")
    # else:
    #    log.info("not enough flats yet need 6")
