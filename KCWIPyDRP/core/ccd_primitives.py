from KCWIPyDRP import PrimitivesBASE
import numpy as np
import pylab as pl


class CcdPrimitives(PrimitivesBASE):

    def __init__(self):
        super(CcdPrimitives, self).__init__()

    def subtract_bias(self):
        tab = self.n_proctab(targtype='MBIAS')
        self.log.info("%d master bias frames found" % len(tab))
        self.img_subtract(tab, suffix='mbias', indir='redux')
        self.frame.header['BIASSUB'] = True
        logstr = self.subtract_bias.__module__ + "." + \
                 self.subtract_bias.__qualname__
        self.frame.header['HISTORY'] = logstr
        self.log.info(logstr)

    def subtract_oscan(self):
        bsec, dsec, tsec, direc = self.map_ccd()
        namps = len(bsec)
        # parameters
        oscanbuf = 20
        minoscan = 75
        if namps == 4:
            porder = 2
        else:
            porder = 7
        # loop over amps
        for ia in range(namps):
            if (bsec[ia][3] - bsec[ia][2]) > minoscan:
                x0 = bsec[ia][2] + oscanbuf
                x1 = bsec[ia][3] - oscanbuf
                y0 = bsec[ia][0]
                y1 = bsec[ia][1] + 1
                osvec = np.nanmedian(self.frame.data[y0:y1, x0:x1], axis=1)
                xx = np.arange(len(osvec), dtype=np.float)
                # fit it
                if direc[ia]:
                    # forward read skips first 50 px
                    oscoef = np.polyfit(xx[50:], osvec[50:], porder)
                else:
                    # reverse read skips last 50 px
                    oscoef = np.polyfit(xx[:-50], osvec[:-50], porder)
                # plot it
                osfit = np.polyval(oscoef, xx)
                pl.plot(osvec)
                legend = ["oscan", ]
                pl.plot(osfit)
                legend.append("fit")
                pl.xlabel("pixel")
                pl.ylabel("DN")
                pl.legend(legend)
                pl.title("img #%d amp #%d" % (self.frame.header['FRAMENO'],
                                              (ia+1)))
                pl.show()
                # subtract it
                for ix in range(dsec[ia][2], dsec[ia][3]):
                    self.frame.data[y0:y1, ix] -= osfit
                self.frame.header['OSCANSUB'] = True
            else:
                self.log.info("not enough overscan to fit amp %d")
                self.frame.header['OSCANSUB'] = False

        logstr = self.subtract_oscan.__module__ + "." + \
                 self.subtract_oscan.__qualname__
        self.frame.header['HISTORY'] = logstr
        self.log.info(logstr)

    def trim_oscan(self):
        bsec, dsec, tsec, direc = self.map_ccd()
        max_sec = max(tsec)
        new = np.zeros((max_sec[1]+1, max_sec[3]+1), dtype=np.uint16)
        namps = len(dsec)
        for ia in range(namps):
            yi0 = dsec[ia][0]
            yi1 = dsec[ia][1] + 1
            xi0 = dsec[ia][2]
            xi1 = dsec[ia][3] + 1
            yo0 = tsec[ia][0]
            yo1 = tsec[ia][1] + 1
            xo0 = tsec[ia][2]
            xo1 = tsec[ia][3] + 1
            new[yo0:yo1, xo0:xo1] = self.frame.data[yi0:yi1, xi0:xi1]
            sec = "[%d:" % (xo0+1)
            sec += "%d," % xo1
            sec += "%d:" % (yo0+1)
            sec += "%d]" % yo1
            self.frame.header['ASEC%d' % (ia+1)] = sec
            # self.frame.header.pop('ASEC%d' % (ia + 1))
            self.frame.header.pop('BSEC%d' % (ia + 1))
            self.frame.header.pop('DSEC%d' % (ia + 1))
            self.frame.header.pop('CSEC%d' % (ia + 1))

        self.frame.data = new
        self.frame.header['NAXIS1'] = max_sec[3] + 1
        self.frame.header['NAXIS2'] = max_sec[1] + 1

        logstr = self.trim_oscan.__module__ + "." + \
                 self.trim_oscan.__qualname__
        self.frame.header['HISTORY'] = logstr
        self.log.info(logstr)

    def correct_gain(self):
        self.log.info("correct_gain")

    def remove_crs(self):
        self.log.info("remove_crs")

    def remove_badcols(self):
        self.log.info("remove_badcols")

    def rectify_image(self):
        self.log.info("rectify_image")

    def parse_imsec(self, section_key=None):
        # forward read?
        xfor = True
        yfor = True
        section = self.frame.header[section_key]
        p1 = int(section[1:-1].split(',')[0].split(':')[0])
        p2 = int(section[1:-1].split(',')[0].split(':')[1])
        p3 = int(section[1:-1].split(',')[1].split(':')[0])
        p4 = int(section[1:-1].split(',')[1].split(':')[1])
        # tests for individual axes
        if p1 > p2:
            x0 = p2 - 1
            x1 = p1 - 1
            xfor = False
        else:
            x0 = p1 - 1
            x1 = p2 - 1
        if p3 > p4:
            y0 = p4 - 1
            y1 = p3 - 1
            yfor = False
        else:
            y0 = p3 - 1
            y1 = p4 - 1
        # package output
        sec = (y0, y1, x0, x1)
        rfor = (yfor, xfor)
        # use python axis ordering
        return sec, rfor

    def map_ccd(self):
        """Return CCD section variables useful for processing

        Uses FITS keyword NVIDINP to determine how many amplifiers were used
        to read out the CCD.  Then reads the corresponding BSECn, and
        DSECn keywords, where n is the amplifier number.  The indices are
        converted to Python (0-biased, y axis first) indices and an array
        is constructed for each of the two useful sections of the CCD as
        follows:

        Bsec[0][0] - First amp, y lower limit
        Bsec[0][1] - First amp, y upper limit
        Bsec[0][2] - First amp, x lower limit
        Bsec[0][3] - First amp, x upper limit
        Bsec[1][0] - Second amp, y lower limit
        etc.

        Bsec is the full overscan region for the given amplifier and is used
        to calculate and perform the overscan subtraction.

        Dsec is the full CCD region for the given amplifier and is used to
        trim the image after overscan subtraction has been performed.

        Tsec accounts for trimming the image according to Dsec.

        Amps are assumed to be organized as follows:

        (0,ny)	--------- (nx,ny)
                | 3 | 4 |
                ---------
                | 1 | 2 |
        (0,0)	--------- (nx, 0)

        Args:
            self: instance of CcdPrimitive class (automatic)

        Returns:
            list: (int) y0, y1, x0, x1 for bias section
            list: (int) y0, y1, x0, x1 for data section
            list: (int) y0, y1, x0, x1 for trimmed section
            list: (bool) y-direction, x-direction, True if forward, else False

        """
        namps = self.frame.header['NVIDINP']
        # TODO: check namps
        # section lists
        bsec = []
        dsec = []
        tsec = []
        direc = []
        # loop over amps
        for i in range(namps):
            sec, rfor = self.parse_imsec(section_key='BSEC%d' % (i+1))
            bsec.append(sec)
            sec, rfor = self.parse_imsec(section_key='DSEC%d' % (i+1))
            dsec.append(sec)
            direc.append(rfor)
            if i == 0:
                y0 = 0
                y1 = sec[1] - sec[0]
                x0 = 0
                x1 = sec[3] - sec[2]
            elif i == 1:
                y0 = 0
                y1 = sec[1] - sec[0]
                x0 = tsec[0][3] + 1
                x1 = x0 + sec[3] - sec[2]
            elif i == 2:
                y0 = tsec[0][1] + 1
                y1 = y0 + sec[1] - sec[0]
                x0 = 0
                x1 = sec[3] - sec[2]
            elif i == 3:
                y0 = tsec[0][1] + 1
                y1 = y0 + sec[1] - sec[0]
                x0 = tsec[0][3] + 1
                x1 = x0 + sec[3] - sec[2]
            else:
                # should not get here
                y0 = -1
                y1 = -1
                x0 = -1
                x1 = -1
                self.log.info("ERROR - bad amp number: %d" % i)
            tsec.append((y0, y1, x0, x1))

        return bsec, dsec, tsec, direc

