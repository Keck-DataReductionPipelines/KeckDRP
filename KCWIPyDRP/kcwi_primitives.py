from KCWIPyDRP.ccd_primitives import CcdPrimitives
from KCWIPyDRP.imgmath_primitives import ImgmathPrimitives
from astropy import log
from astropy.table import Table
import os


class KcwiPrimitives(CcdPrimitives, ImgmathPrimitives):

    def set_frame(self, frame):
        self.frame = frame

    def generate_output_image(self, suffix=None, outdir='redux'):
        if suffix is not None:
            origfn = self.frame.header['OFNAME']
            outfn = os.path.join(outdir,
                                 origfn.split('.')[0]+'_'+suffix+'.fits')
            self.frame.write(outfn)
            log.info("output file: %s" % outfn)

    def subtract_scattered_light(self):
        log.info("subtract_scattered_light")

    def new_proctab(self):
        cnames = ('CID', 'TYPE', 'CAM', 'GRAT', 'GANG', 'CWAVE', 'BIN', 'FILT',
                  'MJD', 'FRAMENO', 'STAGE', 'SUFF', 'OFNAME')
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

    def update_proctab(self, suffix='raw', newtype=None):
        if self.frame is not None and self.proctab is not None:
            stages = {'raw': 0,
                      'int': 1,
                      'intd': 2,
                      'intf': 3,
                      'intk': 4,
                      'cube': 5,
                      'cubed': 6,
                      'cubes': 7}
            if suffix in stages:
                stage = stages[suffix]
            else:
                stage = 9
            if newtype is not None:
                self.frame.header['IMTYPE'] = newtype
            # output file
            # if stage > 0:
            #    self.generate_output_image(suffix=suffix)
            # new row for proc table
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
                       stage,
                       suffix,
                       self.frame.header['OFNAME']]
        else:
            new_row = None
        self.proctab.add_row(new_row)

    def n_proctab(self, targtype=None):
        if targtype is not None and self.proctab is not None:
            tab = self.proctab[(self.proctab['TYPE'] == targtype)]
            tab = tab[(tab['CID'] == self.frame.header['STATEID'])]
        else:
            tab = None

        return tab
