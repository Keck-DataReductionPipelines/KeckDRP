from astropy.nddata import CCDData, NDData
from .KCWI import KcwiConf


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

    def camang(self):
        if self.camera() == 0:      # Blue
            key = 'BARTANG'
        elif self.camera() == 1:    # Red
            key = 'RARTANG'
        else:
            raise ValueError("unable to determine camera angle: "
                             "CAMERA undefined")
        return self.header[key]

    def filter(self):
        if self.camera() == 0:      # Blue
            filt = self.header['BFILTNAM']
        elif self.camera() == 1:    # Red
            filt = 'None'
        else:
            raise ValueError("unable to determine filter: "
                             "CAMERA undefined")
        return filt

    def grangle(self):
        if self.camera() == 0:      # Blue
            key = 'BGRANGLE'
        elif self.camera() == 1:    # Red
            key = 'RGRANGLE'
        else:
            raise ValueError("unable to determine grating angle: "
                             "CAMERA undefined")
        return self.header[key]

    def grating(self):
        if self.camera() == 0:      # Blue
            key = 'BGRATNAM'
        elif self.camera() == 1:    # Red
            key = 'RGRATNAM'
        else:
            raise ValueError("unable to determine grating: CAMERA undefined")
        return self.header[key]

    def adjang(self):
        if 'BH' in self.grating() or 'RH' in self.grating():
            return 180.0
        if 'BM' in self.grating() or 'RM' in self.grating():
            return 0.0
        if 'BL' in self.grating() or 'RL' in self.grating():
            return 0.0

    def rho(self):
        if 'BH1' in self.grating():
            return 3.751
        elif 'BH2' in self.grating():
            return 3.255
        elif 'BH3' in self.grating():
            return 2.800
        elif 'RH1' in self.grating():
            return 2.420
        elif 'RH2' in self.grating():
            return 2.030
        elif 'RH3' in self.grating():
            return 1.705
        elif 'RH4' in self.grating():
            return 1.435
        elif 'BM' in self.grating():
            return 1.900
        elif 'RM1' in self.grating():
            return 1.220
        elif 'RM2' in self.grating():
            return 0.921
        elif 'BL' in self.grating():
            return 0.870
        elif 'RL' in self.grating():
            return 0.514
        else:
            raise ValueError("unable to compute rho: grating undefined")

    def cwave(self):
        if self.camera() == 0:      # Blue
            key = 'BCWAVE'
        elif self.camera() == 1:    # Red
            key = 'RCWAVE'
        else:
            raise ValueError("unable to determine central wavelength: "
                             "CAMERA undefined")
        return self.header[key]

    def resolution(self, refwave=None):
        """Return FWHM resolution in Angstroms for the given grating"""
        # get reference wavelength
        if refwave:
            rw = refwave
        else:
            if 'B' in self.grating():
                rw = 4500.
            else:
                rw = 7500.
        # Calc rez from grating resolution (d-lambda = lambda/R)
        # First, assume large slicer IFU
        if 'BH' in self.grating():
            rez = rw / 5000.
        elif 'RH' in self.grating():
            rez = rw / 5000.
        elif 'BM' in self.grating():
            rez = rw / 2500.
        elif 'RM' in self.grating():
            rez = rw / 2500.
        elif 'BL' in self.grating():
            rez = rw / 1250.
        elif 'RL' in self.grating():
            rez = rw / 1250.
        else:
            raise ValueError("unable to compute atlas resolution: "
                             "grating undefined")
        # Adjust for slicer
        if self.ifunum() == 2:      # Medium slicer
            rez /= 2.
        elif self.ifunum() == 3:    # Small slicer
            rez /= 4.

        return rez

    def namps(self):
        return self.header['NVIDINP']

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

    def minoscanpix(self):
        return KcwiConf.MINOSCANPIX

    def oscanbuf(self):
        return KcwiConf.OSCANBUF

    def plotlabel(self):
        lab = "Img # %d " % self.header['FRAMENO']
        lab += "(%s) " % self.illum()
        lab += "Sl: %s " % self.ifuname()
        lab += "Fl: %s " % self.filter()
        lab += "Gr: %s " % self.grating()
        return lab

    def plotpause(self):
        return KcwiConf.PLOTPAUSE

    def saveintims(self):
        return KcwiConf.SAVEINTIMS

    def inter(self):
        return KcwiConf.INTER

    def ifuname(self):
        return self.header['IFUNAM']

    def ifunum(self):
        return self.header['IFUNUM']

    def imtype(self):
        return self.header['IMTYPE']

    def illum(self):
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
        return self.header['ILLUM']

    def taperfrac(self):
        return KcwiConf.TAPERFRAC