from atomict.datasets import get


import pytest


def test_get_success():

    df = get('clintox')
    assert df.shape == (1484, 3)


def test_get_failure():

    with pytest.raises(ValueError):
        get('xxx')
