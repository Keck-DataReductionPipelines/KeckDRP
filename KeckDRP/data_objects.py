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
