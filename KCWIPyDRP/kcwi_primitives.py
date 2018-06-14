from KCWIPyDRP.ccd_primitives import CcdPrimitives
from KCWIPyDRP.imgmath_primitives import ImgmathPrimitives
from astropy import log
from astropy.table import Table
import os


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives):

    def set_frame(self, frame):
        self.frame = frame

    def output_master(self, master_type="BIAS"):
        log.info("output_master %s" % master_type)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")

    def new_proctab(self):
        cnames = ('CID', 'TYPE', 'CAM', 'GRAT', 'GANG', 'CWAVE', 'BIN', 'FILT',
                  'MJD', 'FRAMENO', 'SUFF', 'STAGE', 'OFNAME')
        dtypes = ('S24', 'S9', 'S4', 'S5', 'float64', 'float64', 'S4', 'S5',
                  'float64', 'int32', 'int32', 'S5', 'S25')
        meta = {'KCWI DRP PROC TABLE': 'new table'}
        self.proctab = Table(names=cnames, dtype=dtypes, meta=meta)
        # prevent string column truncation
        for col in self.proctab.itercols():
            if col.dtype.kind in 'SU':
                self.proctab.replace_column(col.name, col.astype('object'))

    def read_proctab(self, tfil='kcwi.proc'):
        if os.path.isfile(tfil):
            log.info("reading proc table file: %s" % tfil)
            self.proctab = Table.read(tfil, format='ascii')
        else:
            log.info("proc table file not found: %s" % tfil)
            self.new_proctab()
        # format columns
        self.proctab['GANG'].format = '7.2f'
        self.proctab['CWAVE'].format = '8.2f'
        self.proctab['MJD'].format = '15.6f'
        # prevent string column truncation
        for col in self.proctab.itercols():
            if col.dtype.kind in 'SU':
                self.proctab.replace_column(col.name, col.astype('object'))

    def write_proctab(self, tfil='kcwi.proc'):
        if self.proctab is not None:
            if os.path.isfile(tfil):
                over_write = True
            else:
                over_write = False

            self.proctab.write(filename=tfil, format='ascii',
                               overwrite=over_write)
            log.info("writing proc table file: %s" % tfil)
        else:
            log.info("no proc table to write")

    def row_proctab(self, suffix='int'):
        if self.frame is not None and self.proctab is not None:
            stages = {'int': 1,
                      'intd': 2,
                      'intf': 3,
                      'intk': 4,
                      'cube': 5,
                      'cubed': 6,
                      'cubes': 7}
            new_row = [self.frame.header['STATEID'],
                       self.frame.header['IMTYPE'],
                       self.frame.header['CAMERA'],
                       self.frame.header['BGRATNAM'],
                       self.frame.header['BGRANGLE'],
                       self.frame.header['BCWAVE'],
                       self.frame.header['BINNING'],
                       self.frame.header['BFILTNAM'],
                       self.frame.header['MJD'],
                       self.frame.header['FRAMENO'],
                       stages[suffix],
                       suffix,
                       self.frame.header['OFNAME']]
        else:
            new_row = None
        return new_row

    def updt_proctab(self, row=None):
        if row is not None:
            self.proctab.add_row(row)

    def in_proctab(self, row=None):
        if row is not None and self.proctab is not None:
            tlen = len(self.proctab)
            # do real tests in here
            return False
        else:
            return False



