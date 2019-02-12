from KeckDRP import PrimitivesBASE
import KeckDRP
import os

import ccdproc


class ImgmathPrimitives(PrimitivesBASE):

    def image_combine(self, tab=None, combine_type='bias', in_directory=None,
                      suffix=None, method='average', unit='adu', keylog=None):
        if tab is not None:
            file_list = tab['OFNAME']
            image_numbers = tab['FRAMENO']

            if in_directory is None:
                prefix = '.'
            else:
                prefix = in_directory

            if suffix is None:
                suffix = '.fits'
            else:
                suffix = '_' + suffix + '.fits'

            # stack images
            stack = []
            for f in file_list:
                infile = os.path.join(prefix, f.split('.')[0] + suffix)
                self.log.info("reading image: %s" % infile)
                stack.append(KeckDRP.KcwiCCD.read(infile, unit=unit))
            # combine biases
            if 'bias' in combine_type:
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
            # handle missing CCDCFG
            if 'CCDCFG' not in self.frame.header:
                ccdcfg = self.frame.header['CCDSUM'].replace(" ", "")
                ccdcfg += "%1d" % self.frame.header['CCDMODE']
                ccdcfg += "%02d" % self.frame.header['GAINMUL']
                ccdcfg += "%02d" % self.frame.header['AMPMNUM']
                self.frame.header['CCDCFG'] = ccdcfg
            # do we log list in header?
            if keylog is not None:
                # make a list of image numbers string
                image_numbers_string = ','.join(str(e) for e in image_numbers)
                card = (image_numbers_string, '')
                # was a keyword comment provided?
                if keylog is not None:
                    if keylog in self.keyword_comments:
                        card = (image_numbers_string,
                                self.keyword_comments[keylog])
                    else:
                        card = (image_numbers_string, '')
                self.frame.header[keylog] = card
            log_string = self.image_combine.__module__ + "." + \
                         self.image_combine.__qualname__
            self.frame.header['HISTORY'] = log_string
            self.log.info("%s %s using %s" % (self.image_combine.__name__,
                                              combine_type, method))
        else:
            self.log.error("something went wrong with image_combine")

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
                subtrahend = KeckDRP.KcwiCCD.read(infile, unit=unit)

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
                self.log.error("Error - empty subtrahend")
        else:
            self.log.error("no image found to subtract")

    def img_divide(self):
        self.log.info("img_divide")
