import json
import warnings

from .msgpack import atoms_to_dict as msgpack_atoms_to_dict, dict_to_atoms as msgpack_dict_to_atoms

try:
    from ase import Atoms
    from ase.io import read, write
except ImportError:
    raise ImportError(
        "The 'ase' package is required for JSON operations with Atoms objects. "
        "To install the optional dependencies such as ase, spglib, pymatgen, use `pip install atomict[utils]`"
    )


DEP_WARNING = """
This function is deprecated. Use atomict.io.msgpack.atoms_to_dict instead.
"""

def atoms_to_json(atoms: Atoms) -> str:
    warnings.warn(DEP_WARNING, DeprecationWarning, stacklevel=2)
    encoded = msgpack_atoms_to_dict([atoms])
    return json.dumps(encoded[0])


def atoms_to_dict(atoms: Atoms) -> dict:
    warnings.warn(DEP_WARNING, DeprecationWarning, stacklevel=2)
    encoded = msgpack_atoms_to_dict([atoms])
    return encoded[0]


def json_to_atoms(json_str: str) -> Atoms:
    warnings.warn(DEP_WARNING, DeprecationWarning, stacklevel=2)
    encoded = json.loads(json_str)
    return msgpack_dict_to_atoms([encoded])[0]
