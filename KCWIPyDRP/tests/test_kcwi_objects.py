from astropy.nddata import CCDData
from KCWIPyDRP.kcwi_objects import KcwiCCD
import pytest
import numpy as np
from astropy.utils import NumpyRNGContext
from astropy import units as u

DEFAULTS = {
    'seed': 123,
    'data_size': 100,
    'data_scale': 1.0,
    'data_mean': 0.0
}

DEFAULT_SEED = 123
DEFAULT_DATA_SIZE = 100
DEFAULT_DATA_SCALE = 1.0


def value_from_markers(key, request):
    try:
        val = request.keywords[key].args[0]
    except KeyError:
        val = DEFAULTS[key]
    return val


@pytest.fixture
def ccd_data(request):
    size = value_from_markers('data_size', request)
    scale = value_from_markers('data_scale', request)
    mean = value_from_markers('data_mean', request)

    with NumpyRNGContext(DEFAULTS['seed']):
        data = np.random.normal(loc=mean, size=[size, size], scale=scale)

    fake_meta = {'my_key': 42, 'your_key': 'not 42'}
    ccd = CCDData(data, unit=u.adu)
    ccd.header = fake_meta
    return ccd


def test_ccddata_empty():
    with pytest.raises(TypeError):
        KcwiCCD()  # empty initializer should fail


def test_ccddata_must_have_unit():
    with pytest.raises(ValueError):
        KcwiCCD(np.zeros([100, 100]))


def test_ccddata_unit_cannot_be_set_to_none(ccd_data):
    with pytest.raises(TypeError):
        ccd_data.unit = None


def test_ccddata_meta_header_conflict():
    with pytest.raises(ValueError) as exc:
        KcwiCCD([1, 2, 3], unit='', meta={1: 1}, header={2: 2})
        assert "can't have both header and meta." in str(exc)
