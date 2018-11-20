# from KCWIPyDRP.ccd_primitives import CcdPrimitives
# from KCWIPyDRP.imgmath_primitives import ImgmathPrimitives
from .. import kcwi_primitives
from .. import kcwi_objects
# from astropy import log
# from astropy.table import Table
import numpy as np
# import os
import pytest


@pytest.fixture
def p():
    p = kcwi_primitives.KcwiPrimitives()
    return p


def test_set_frame(p):
    myframe = kcwi_objects.KcwiCCD(np.random.normal(size=(10, 10)), unit="adu")
    p.set_frame(myframe)
    assert p.frame == myframe
