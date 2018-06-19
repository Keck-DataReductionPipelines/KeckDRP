from KCWIPyDRP import PrimitivesBASE
from ..core import CcdPrimitives
from ..core import ImgmathPrimitives
from ..core import ProctabPrimitives
import os


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives, ProctabPrimitives):

    def __init__(self):
        super(KcwiPrimitives, self).__init__()

    def write_image(self, suffix=None, outdir='redux'):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(outdir,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            self.frame.write(outfn)
            self.log.info("output file: %s" % outfn)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")
