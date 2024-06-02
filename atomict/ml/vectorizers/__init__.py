import numpy as np
from rdkit import Chem

from atomict.ml.vectorizers.bitvect import avalon_fingerprints, morgan_fingerprints
from atomict.ml.vectorizers.domainrelevant import domainrelevant_fingerprint


def get_amd_fingerprint(smiles: str, remove_extraneous_hydrogens: bool = False, morgan_radius: int = 2):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    if remove_extraneous_hydrogens:
        mol = Chem.RemoveHs(mol)

    avalon = avalon_fingerprints(mol)
    morgan = morgan_fingerprints(mol, radius=morgan_radius)
    domain_r = domainrelevant_fingerprint(mol)

    return np.concatenate([avalon, morgan, domain_r], axis=0)
