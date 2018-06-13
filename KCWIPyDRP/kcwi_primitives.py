from KCWIPyDRP.ccd import Ccd
from KCWIPyDRP.imgmath import Imgmath
from astropy import log
from astropy.nddata import CCDData


class Kcwi(Ccd, Imgmath):

    img = None

    def __init__(self, imgf):
        Ccd.__init__(self)
        Imgmath.__init__(self)
        # read image
        self.img = CCDData.read(imgf, unit="adu")
        # set ILLUM keyword
        if self.img.header['IMTYPE'] == 'ARCLAMP':
            if self.img.header['LMP0STAT'] == 1 and self.img.header['LMP0SHST'] == 1:
                self.img.header['ILLUM'] = self.img.header['LMP0NAM']
            elif self.img.header['LMP1STAT'] == 1 and self.img.header['LMP1SHST'] == 1:
                self.img.header['ILLUM'] = self.img.header['LMP1NAM']
            else:
                self.img.header['ILLUM'] = 'Test'
        elif self.img.header['IMTYPE'] == 'FLATLAMP':
            if self.img.header['LMP3STAT'] == 1:
                self.img.header['ILLUM'] = 'Contin'
            else:
                self.img.header['ILLUM'] = 'Test'
        elif self.img.header['IMTYPE'] == 'CONTBARS':
            if self.img.header['LMP3STAT'] == 1:
                self.img.header['ILLUM'] = 'Contin'
            else:
                self.img.header['ILLUM'] = 'Test'
        elif self.img.header['IMTYPE'] == 'OBJECT':
            self.img.header['ILLUM'] = 'Object'
        else:
            self.img.header['ILLUM'] = 'Test'
        # log results
        log.info(self.img.header['ILLUM'])
        log.info("init")

    def output_master(self, master_type="BIAS"):
        log.info("output_master %s" % master_type)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")
