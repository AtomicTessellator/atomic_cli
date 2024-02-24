import numpy as np
from rdkit import Chem

from .bitvect import avalon_fingerprints, morgan_fingerprints
from .domainrelevant import domainrelevant_fingerprint


# Better name, get_full_fingerprint will be deprecated
def get_amd_fingerprint(smiles: str, remove_extraneous_hydrogens: bool = False):
    return get_full_fingerprint(smiles, remove_extraneous_hydrogens)


def get_full_fingerprint(smiles: str, remove_extraneous_hydrogens: bool = False):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    if remove_extraneous_hydrogens:
        mol = Chem.RemoveHs(mol)

    avalon = avalon_fingerprints(mol)
    morgan = morgan_fingerprints(mol)
    domain_r = domainrelevant_fingerprint(mol)

    return np.concatenate([avalon, morgan, domain_r], axis=0)
