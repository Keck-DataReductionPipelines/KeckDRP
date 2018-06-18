from astropy import log


class CcdPrimitives:

    def set_frame(self, frame):
        self.frame = frame

    def subtract_bias(self):
        tab = self.n_proctab(targtype='MBIAS')
        log.info("%d master bias frames found" % len(tab))
        self.img_subtract(tab, suffix='mbias', indir='redux')

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

    def parse_imsec(self, section_key=None):
        section = self.frame.header[section_key]
        p1 = int(section[1:-1].split(',')[0].split(':')[0])
        p2 = int(section[1:-1].split(',')[0].split(':')[1])
        p3 = int(section[1:-1].split(',')[1].split(':')[0])
        p4 = int(section[1:-1].split(',')[1].split(':')[1])
        # tests
        if p1 > p2:
            x0 = p2
            x1 = p1
        else:
            x0 = p1
            x1 = p2
        if p3 > p4:
            y0 = p4
            y1 = p3
        else:
            y0 = p3
            y1 = p4
        # use python axis ordering
        return y0, y1, x0, x1
