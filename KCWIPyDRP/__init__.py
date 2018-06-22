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
from .lookups import keyword_comments

# set up namespace, unless we are in setup...
if not _ASTROPY_SETUP_:
    # from .core import *
    # from .ccddata import *
    # from .combiner import *
    # from .image_collection import *
    from astropy import config as _config

    class Conf(_config.ConfigNamespace):
        """
        Configuration parameters for KCWIPyDRP.
        """
        CRZAP = _config.ConfigItem(
            True,
            'Perform cosmic ray rejection'
            )
        INTER = _config.ConfigItem(
            True,
            'Interactive operation'
            )
        PLOTPAUSE = _config.ConfigItem(
            0.5,
            'Pause length between plots in seconds'
            )
        MINOSCANPIX = _config.ConfigItem(
            75,
            'Minimum number of pixels for overscan'
            )
        OSCANBUF = _config.ConfigItem(
            20,
            'Pixel buffer to exclude at edges of overscan'
            )
    conf = Conf()


class PrimitivesBASE():
    def __init__(self):
        self.frame = None
        self.log = log
        self.log.enable_color()
        self.conf = conf
        self.keyword_comments = keyword_comments.keyword_comments

    def set_frame(self, frame):
        self.frame = frame
