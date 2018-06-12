from KCWIPyDRP.ccd import ccd
from KCWIPyDRP.imgmath import imgmath
from astropy import log

class kcwi(ccd,imgmath):
    def __init__(self):
        super(kcwi,self).__init__()
        log.info("init")
    def output_master(self,master_type="BIAS"):
        log.info("output_master %s" % master_type)
    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")

