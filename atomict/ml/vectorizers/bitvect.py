import numpy as np
from rdkit import DataStructs
from rdkit.Avalon.pyAvalonTools import GetAvalonCountFP
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Mol


def _to_np(fp):
    d_np = np.zeros(0)
    DataStructs.ConvertToNumpyArray(fp, d_np)
    return d_np


def avalon_fingerprints(mol: Mol, n_bits: int = 2048):
    descript = GetAvalonCountFP(mol, nBits=n_bits)
    return _to_np(descript)


def morgan_fingerprints(mol: Mol, n_bits: int = 2048, radius: int = 2):
    descript = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return _to_np(descript)
