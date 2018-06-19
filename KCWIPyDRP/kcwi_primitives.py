from astropy.table import Table


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives):

    def __init__(self):
        super(KcwiPrimitives, self).__init__()

    def write_image(self, suffix=None, outdir='redux'):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(outdir,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            self.frame.write(outfn)
            log.info("output file: %s" % outfn)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")
