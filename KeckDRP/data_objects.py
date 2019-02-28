from astropy.nddata import CCDData


class KcwiCCD(CCDData):
    """
    the KCWICCD class subclasses the CCDData class, which is subclass of NDData
    """
    def __init__(self, *args, **kwd):
        if 'meta' not in kwd:
            kwd['meta'] = kwd.pop('header', None)
        if 'header' in kwd:
            raise ValueError("can't have both header and meta.")
        super().__init__(*args, **kwd)

    def camera(self):
        if 'CAMERA' in self.header:
            if 'BLUE' in self.header['CAMERA']:
                return 0
            elif 'RED' in self.header['CAMERA']:
                return 1
            else:
                return -1
        else:
            return -1

    def nasmask(self):
        if self.camera() == 0:      # Blue
            if 'Mask' in self.header['BNASNAM']:
                return True
            else:
                return False
        elif self.camera() == 1:    # Red
            if 'Mask' in self.header['RNASNAM']:
                return True
            else:
                return False
        else:
            raise ValueError("unable to determine mask: CAMERA undefined")

    def xbinsize(self):
        return int(self.header['BINNING'].split(',')[0])

    def ybinsize(self):
        return int(self.header['BINNING'].split(',')[-1])

    def imtype(self):
        # set ILLUM keyword
        # ARCS
        if self.header['IMTYPE'] == 'ARCLAMP':
            if self.header['LMP0STAT'] == 1 and \
                    self.header['LMP0SHST'] == 1:
                self.header['ILLUM'] = self.header['LMP0NAM']
            elif self.header['LMP1STAT'] == 1 and \
                    self.header['LMP1SHST'] == 1:
                self.header['ILLUM'] = self.header['LMP1NAM']
            else:
                self.header['ILLUM'] = 'Test'
        # Internal FLATS
        elif self.header['IMTYPE'] == 'FLATLAMP':
            if self.header['LMP3STAT'] == 1:
                self.header['ILLUM'] = 'Contin'
            else:
                self.header['ILLUM'] = 'Test'
        # DOMES
        elif self.header['IMTYPE'] == 'DOMEFLAT':
            if self.header['FLIMAGIN'] == 'on' or \
                    self.header['FLSPECTR'] == 'on':
                self.header['ILLUM'] = 'Dome'
            else:
                self.header['ILLUM'] = 'Test'
        # Twilight FLATS
        elif self.header['IMTYPE'] == 'TWIFLAT':
            self.header['ILLUM'] = 'Twilit'
        # BARS
        elif self.header['IMTYPE'] == 'CONTBARS':
            if self.header['LMP3STAT'] == 1:
                self.header['ILLUM'] = 'Contin'
            else:
                self.header['ILLUM'] = 'Test'
        # OBJECT
        elif self.header['IMTYPE'] == 'OBJECT':
            self.header['ILLUM'] = 'Object'
        else:
            self.header['ILLUM'] = 'Test'
