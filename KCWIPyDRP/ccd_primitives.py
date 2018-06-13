import numpy as np
import scipy as sp
from astropy import log


class CcdPrimitives:

    def set_frame(self,frame):
        self.frame = frame

    def subtract_oscan(self):
        log.info(self.frame.header['BSEC1'])
        log.info("subtract_oscan")

    def trim_oscan(self):
        log.info("trim_oscan")

    def correct_gain(self):
        log.info("correct_gain")

    def remove_crs(self):
        log.info("remove_crs")

    def remove_badcols(self):
        log.info("remove_badcols")

    def rectify_image(self):
        log.info("rectify_image")

