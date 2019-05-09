# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This packages contains the KPF specific primitives and recipes
"""
from . import recipes

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
        False,
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
    MINIMUM_NUMBER_OF_BIASES = _config.ConfigItem(
        7,
        'Minimum number of biases'
    )
    MINIMUM_NUMBER_OF_DARKS = _config.ConfigItem(
        3,
        'Minimum number of darks'
    )
    MINIMUM_NUMBER_OF_FLATS = _config.ConfigItem(
        6,
        'Minimum number of flats'
    )


KpfConf = Conf()
