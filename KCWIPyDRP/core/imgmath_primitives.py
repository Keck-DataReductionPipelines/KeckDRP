from KCWIPyDRP import PrimitivesBASE
import os
from ..kcwi.kcwi_objects import KcwiCCD
import ccdproc


class ImgmathPrimitives(PrimitivesBASE):

    def img_combine(self, tab=None, ctype='bias', indir=None, suffix=None,
                    method='average', unit='adu', keylog=None):
        if tab is not None:
            flist = tab['OFNAME']
            imnos = tab['FRAMENO']

            if indir is None:
                pref = '.'
            else:
                pref = indir

            if suffix is None:
                suff = '.fits'
            else:
                suff = '_' + suffix + '.fits'

            # stack images
            stack = []
            for f in flist:
                infile = os.path.join(pref, f.split('.')[0] + suff)
                self.log.info("reading image: %s" % infile)
                stack.append(KcwiCCD.read(infile, unit=unit))
            # combine biases
            if 'bias' in ctype:
                self.set_frame(ccdproc.combine(stack, method=method,
                                               sigma_clip=True,
                                               sigma_clip_low_thresh=None,
                                               sigma_clip_high_thresh=2.0))
            # or combine any other type
            else:
                self.set_frame(ccdproc.combine(stack, method=method))
            self.frame.header['NSTACK'] = (len(stack),
                                           self.keyword_comments['NSTACK'])
            self.frame.header['STCKMETH'] = (method,
                                             self.keyword_comments['STCKMETH'])
            # do we log list in header?
            if keylog is not None:
                # make a list of image numbers string
                imnos_str = ','.join(str(e) for e in imnos)
                card = (imnos_str, '')
                # was a keyword comment provided?
                if keylog is not None:
                    if keylog in self.keyword_comments:
                        card = (imnos_str, self.keyword_comments[keylog])
                    else:
                        card = (imnos_str, '')
                self.frame.header[keylog] = card
            logstr = self.img_combine.__module__ + "." + \
                     self.img_combine.__qualname__
            self.frame.header['HISTORY'] = logstr
            self.log.info("%s %s using %s" % (self.img_combine.__name__,
                                              ctype, method))
        else:
            self.log.info("something went wrong with img_combine")

    def img_subtract(self, tab=None, indir=None, suffix=None, unit='adu',
                     keylog=None):
        if tab is not None:
            flist = tab['OFNAME']

            if indir is None:
                pref = '.'
            else:
                pref = indir

            if suffix is None:
                suff = '.fits'
            else:
                suff = '_' + suffix + '.fits'

            subtrahend = None
            for f in flist:
                infile = os.path.join(pref, f.split('.')[0] + suff)
                self.log.info("reading image to subtract: %s" % infile)
                subtrahend = KcwiCCD.read(infile, unit=unit)

            if subtrahend is not None:
                result = self.frame.subtract(subtrahend)
                result.meta = self.frame.meta
                self.set_frame(result)
                if keylog is not None:
                    if keylog in self.keyword_comments:
                        card = (infile, self.keyword_comments[keylog])
                    else:
                        card = (infile, '')
                    self.frame.header[keylog] = card
                logstr = self.img_subtract.__module__ + "." + \
                         self.img_subtract.__qualname__
                self.frame.header['HISTORY'] = logstr
                self.log.info(self.img_subtract.__qualname__)
            else:
                self.log.info("Error - empty subtrahend")
        else:
            self.log.info("no image found to subtract")

    def img_divide(self):
        self.log.info("img_divide")
