from KCWIPyDRP.ccd import Ccd
from KCWIPyDRP.imgmath import Imgmath
from astropy import log


class Kcwi(Ccd, Imgmath):

    def __init__(self):
        Ccd.__init__(self)
        Imgmath.__init__(self)
        log.info("init")

    def output_master(self, master_type="BIAS"):
        log.info("output_master %s" % master_type)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")
