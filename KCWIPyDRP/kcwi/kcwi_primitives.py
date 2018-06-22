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

    def subtract_bias(self):
        tab = self.n_proctab(targtype='MBIAS')
        self.log.info("%d master bias frames found" % len(tab))
        self.img_subtract(tab, suffix='mbias', indir='redux', keylog='MBFILE')
        self.frame.header['BIASSUB'] = (True, self.keyword_comments['BIASSUB'])
        logstr = self.subtract_bias.__module__ + "." + \
                 self.subtract_bias.__qualname__
        self.frame.header['HISTORY'] = logstr
        self.log.info(self.subtract_bias.__qualname__)

    def subtract_scattered_light(self):
        self.log.info("subtract_scattered_light")
