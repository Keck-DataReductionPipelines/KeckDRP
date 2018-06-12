from astropy import log


def make_master_flat(p):
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()
    p.img_combine()
    p.output_master(master_type="FLAT")
