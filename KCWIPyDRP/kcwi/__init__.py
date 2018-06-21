# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This packages contains the KCWI specific primitives and recipes
"""
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
        CRZAP= _config.ConfigItem(
            True,
            'Perform cosmic ray rejection'
            )
    conf = Conf()