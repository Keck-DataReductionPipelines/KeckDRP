# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The KCWIPyDRP is the Python version of the KCWI pipeline
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------
from astropy import log




class PrimitivesBASE():
    def __init__(self):
        self.frame = None
        self.log = log
        self.log.enable_color()

    def set_frame(self, frame):
        self.frame = frame
