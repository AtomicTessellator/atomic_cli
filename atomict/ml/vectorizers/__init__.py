import numpy as np
from rdkit import Chem

from .bitvect import avalon_fingerprints, morgan_fingerprints
from .domainrelevant import domainrelevant_fingerprint


def get_full_fingerprint(smiles: str):

    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return None
    
    avalon = avalon_fingerprints(mol)
    morgan = morgan_fingerprints(mol)
    domain_r = domainrelevant_fingerprint(smiles)
    
    return np.concatenate([avalon, morgan, domain_r], axis=1)
