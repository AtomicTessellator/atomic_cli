import io
import json

try:
    from ase import Atoms
    from ase.io import read, write
except ImportError:
    raise ImportError(
        "The 'ase' package is required for JSON operations with Atoms objects. "
        "To install the optional dependencies such as ase, spglib, pymatgen, use `pip install atomict[utils]`"
    )


def atoms_to_json(atoms: Atoms) -> str:
    with io.StringIO() as buf:
        write(buf, atoms, format="json")
        return buf.getvalue()


def atoms_to_dict(atoms: Atoms) -> dict:
    with io.StringIO() as buf:
        write(buf, atoms, format="json")
        return json.loads(buf.getvalue())


def json_to_atoms(json_str: str) -> Atoms:
    with io.StringIO(json_str) as buf:
        return read(buf, format="json")
