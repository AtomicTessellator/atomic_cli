import numpy as np
from rdkit import Chem

from .bitvect import avalon_fingerprints, morgan_fingerprints
from .domainrelevant import domainrelevant_fingerprint


# Better name, get_full_fingerprint will be deprecated
def get_amd_fingerprint(smiles: str, remove_extraneous_hydrogens: bool = False, morgan_radius: int = 2):
    return get_full_fingerprint(smiles, remove_extraneous_hydrogens, morgan_radius)


def get_full_fingerprint(smiles: str, remove_extraneous_hydrogens: bool = False, morgan_radius: int = 2):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    if remove_extraneous_hydrogens:
        mol = Chem.RemoveHs(mol)

    avalon = avalon_fingerprints(mol)
    morgan = morgan_fingerprints(mol, radius=morgan_radius)
    domain_r = domainrelevant_fingerprint(mol)

    return np.concatenate([avalon, morgan, domain_r], axis=0)
