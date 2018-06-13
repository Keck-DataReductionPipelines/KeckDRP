from KCWIPyDRP.ccd_primitives import CcdPrimitives
from KCWIPyDRP.imgmath_primitives import ImgmathPrimitives
from astropy import log
from astropy.nddata import CCDData


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives):

    def set_frame(self,frame):
        self.frame = frame

    def output_master(self, master_type="BIAS"):
        log.info("output_master %s" % master_type)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")
