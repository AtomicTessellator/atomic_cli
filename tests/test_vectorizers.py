import numpy as np
import pytest
from rdkit import Chem

from atomict.ml.vectorizers import (
    avalon_fingerprints,
    morgan_fingerprints,
)


@pytest.fixture
def benzene_mol():
    return Chem.MolFromSmiles("C1=CC=CC=C1")


def test_avalon_fingerprints(benzene_mol):
    fp = avalon_fingerprints(benzene_mol)
    assert isinstance(fp, np.ndarray), "Fingerprint is not a numpy array"
    assert fp.shape == (2048, ), "Fingerprint shape is incorrect"


def test_morgan_fingerprints(benzene_mol):
    fp = morgan_fingerprints(benzene_mol)
    assert isinstance(fp, np.ndarray), "Fingerprint is not a numpy array"
    assert fp.shape == (2048, ), "Fingerprint shape is incorrect"


def test_fingerprint_values():
    """
    Implement tests here to check if the fingerprint values are as expected.
    This requires having a known fingerprint for a given molecule to compare against.
    """
    pass
