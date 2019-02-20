from KeckDRP import PrimitivesBASE
from astropy.table import Table, unique
import os


class ProctabPrimitives(PrimitivesBASE):

    def __init__(self):
        self.proctab = None
        super(ProctabPrimitives, self).__init__()

    def new_proctab(self):
        cnames = ('CID', 'DID', 'TYPE', 'GRPID', 'TTIME', 'CAM', 'GRAT', 'GANG',
                  'CWAVE', 'BIN', 'FILT', 'MJD', 'FRAMENO', 'STAGE', 'SUFF',
                  'OFNAME', 'OBJECT')
        dtypes = ('S24', 'int64', 'S9', 'S12', 'float64', 'S4', 'S5', 'float64',
                  'float64', 'S4', 'S5', 'float64', 'int32', 'int32', 'S5',
                  'S25', 'S25')
        meta = {'KCWI DRP PROC TABLE': 'new table'}
        self.proctab = Table(names=cnames, dtype=dtypes, meta=meta)
        # prevent string column truncation
        for col in self.proctab.itercols():
            if col.dtype.kind in 'SU':
                self.proctab.replace_column(col.name, col.astype('object'))

    def read_proctab(self, tfil='kcwi.proc'):
        if os.path.isfile(tfil):
            self.log.info("reading proc table file: %s" % tfil)
            self.proctab = Table.read(tfil, format='ascii')
            self.proctab.dtypes = ('S24', 'int64', 'S9', 'S12', 'S4', 'S5',
                                   'float64', 'float64', 'S4', 'S5', 'float64',
                                   'int32', 'int32', 'S5', 'S25', 'S25')
        else:
            self.log.info("proc table file not found: %s" % tfil)
            self.new_proctab()
        if len(self.proctab) == 0:
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
            self.log.info("writing proc table file: %s" % tfil)
        else:
            self.log.info("no proc table to write")

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
            # new row for proc table
            if self.frame.header['STATEID'].strip() == '0':
                self.frame.header['STATEID'] = 'NONE'
            if 'GRPID' not in self.frame.header:
                dto = self.frame.header['DATE-OBS']
                fno = self.frame.header['FRAMENO']
                self.frame.header['GRPID'] = "%s-%s" % (dto, fno)
            new_row = [self.frame.header['STATEID'],
                       self.frame.header['CCDCFG'],
                       self.frame.header['IMTYPE'],
                       self.frame.header['GRPID'],
                       self.frame.header['TTIME'],
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
                       self.frame.header['OFNAME'],
                       self.frame.header['OBJECT'].replace(" ", "")]
        else:
            new_row = None
        print("Attempting to add %s" % str(new_row))
        self.proctab.add_row(new_row)
        self.proctab = unique(self.proctab, keys=['CID', 'FRAMENO', 'STAGE'],
                              keep='last')

    def n_proctab(self, target_type=None, target_group=None, nearest=False):
        if target_type is not None and self.proctab is not None:
            self.log.info('Looking for %s frames' % target_type)
            tab = self.proctab[(self.proctab['TYPE'] == target_type)]
            # BIASES must have the same CCDCFG
            if 'BIAS' in target_type:
                self.log.info('Looking for frames with CCDCFG = %s' %
                              self.frame.header['CCDCFG'])
                tab = tab[(tab['DID'] == int(self.frame.header['CCDCFG']))]
                if target_group is not None:
                    tab = tab[(tab['GRPID'] == target_group)]
            # raw DARKS must have the same CCDCFG and TTIME
            elif target_type == 'DARK':
                self.log.info('Looking for frames with CCDCFG = %s and '
                              'TTIME = %f' % (self.frame.header['CCDCFG'],
                                              self.frame.header['TTIME']))
                tab = tab[tab['DID'] == int(self.frame.header['CCDCFG'])]
                tab = tab[tab['TTIME'] == float(self.frame.header['TTIME'])]
                tab = tab[tab['GRPID'] == target_group]
            # MDARKS must have the same CCDCFG, will be scaled to match TTIME
            elif target_type == 'MDARK':
                self.log.info('Looking for frames with CCDCFG = %s' %
                              self.frame.header['CCDCFG'])
                tab = tab[(tab['DID'] == int(self.frame.header['CCDCFG']))]
            else:
                tab = tab[(tab['CID'] == self.frame.header['STATEID'])]
            # Check if nearest entry is requested
            if nearest and len(tab) > 1:
                tfno = self.frame.header['FRAMENO']
                minoff = 99999
                trow = None
                for row in tab:
                    off = abs(row['FRAMENO'] - tfno)
                    if off < minoff:
                        minoff = off
                        trow = row
                if trow is not None:
                    tab = tab[(tab['FRAMENO'] == trow['FRAMENO'])]
        else:
            tab = None
        return tab

    def in_proctab(self):
        imno_list = self.proctab['FRAMENO']
        if self.frame.header['FRAMENO'] in imno_list:
            return True
        else:
            return False
