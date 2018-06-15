import numpy as np
import scipy as sp
from astropy import log
from astropy.table import Table


class ImgmathPrimitives:

    def set_frame(self, frame):
        self.frame = frame

    def img_combine(self, tab=None, suffix=None):
        if tab is not None:
            flist = tab['OFNAME']
            print(flist)
            log.info("img_combine")

    def img_subtract(self):
        log.info("img_subtract")

    def img_divide(self):
        log.info("img_divide")

