from .. import KcwiConf
from KeckDRP import conf



def handle_flats(p, frame):
    # do basic CCD reduction
    p.set_frame(frame)
    p.read_proctab()
    if p.in_proctab():
        p.log.warning("Already processed")
        return
    p.subtract_bias()
    p.subtract_oscan()
    p.trim_oscan()
    p.correct_gain()
    p.remove_badcols()
    p.remove_crs()
    p.rectify_image()

    # update proc table
    p.update_proctab(suffix='int')
    p.write_proctab()

    # output file
    p.write_image(suffix='int')
    p.log.info("flat reduced")

    p.stackInternalFlats()


