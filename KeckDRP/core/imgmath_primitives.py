from KeckDRP import PrimitivesBASE
import KeckDRP
import os

import ccdproc


class ImgmathPrimitives(PrimitivesBASE):

    def add_to_master(self, combine_tab=None, master_tab=None, ctype='bias', indir=None, suffix=None,
                    method='average', unit='adu', keylog=None):
        if combine_tab is not None:
            file_list = combine_tab['OFNAME']
            image_numbers = combine_tab['FRAMENO']

            if indir is None:
                prefix = '.'
            else:
                prefix = indir

            if suffix is None:
                suff = '.fits'
            else:
                suff = '_' + suffix + '.fits'
        n = len(combine_tab)
        if master_tab is not None:

            master_file_list = master_tab['OFNAME']
            master_image_numbers = master_tab['FRAMENO']
            indir = 'redux'
            if indir is None:
                master_prefix = '.'
            else:
                master_prefix = indir
            suffix = 'mbias'
            if suffix is None:
                master_suff = '.fits'
            else:
                master_suff = '_' + suffix + '.fits'
        # Read master
        print(master_file_list)
        print(file_list)
        master_bias = os.path.join(master_prefix, master_file_list[0].split('.')[0] + master_suff)
        # Read new bias
        new_bias = os.path.join(prefix, file_list[-1].split('.')[0]+suff)
        #
        # combine
        master_frame = KeckDRP.KcwiCCD.read(master_bias, unit=unit)
        master_frame = master_frame.multiply((n-1)/n)
        print("Trying to read file: %s" % new_bias)

        new_frame = KeckDRP.KcwiCCD.read(new_bias, unit=unit)
        metadata = new_frame.meta
        new_frame = new_frame.multiply(1/n)
        stack = [new_frame, master_frame]
        # DO THE MATH!!
        result = ccdproc.combine(stack, method=method,
                                           sigma_clip=True,
                                           sigma_clip_low_thresh=None,
                                           sigma_clip_high_thresh=2.0)
        result.meta = metadata
        print(result.header['STATEID'])
        self.set_frame(result)




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
                stack.append(KeckDRP.KcwiCCD.read(infile, unit=unit))
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
            self.log.error("something went wrong with img_combine")

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
                self.log.error("Error - empty subtrahend")
        else:
            self.log.error("no image found to subtract")

    def img_divide(self):
        self.log.info("img_divide")
