from typing import Union, List, Dict, Tuple
from atomict.io.formats.atraj import write_atraj, read_atraj
from atomict.io.formats.tess import write_tess, read_tess


def save_msgpack_trajectory(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str, metadata: Dict = None):
    """Save atoms to a msgpack trajectory file with metadata.
    
    Parameters:
    -----------
    atoms : Atoms or list of Atoms
        The atoms to save
    filename : str
        The output filename
    metadata : dict, optional
        Additional metadata to store with the trajectory
    """

    if filename.endswith('.atraj'):
        write_atraj(atoms, filename, metadata)
    elif filename.endswith('.tess'):
        write_tess(atoms, filename, metadata)
    else:
        extn = filename.split('.')[-1]
        raise ValueError(f"Unsupported file extension: {extn}")


def load_msgpack_trajectory(filename: str) -> Tuple[List['ase.Atoms'], Dict]:
    """Load atoms from a msgpack trajectory file with metadata.
    
    Parameters:
    -----------
    filename : str
        The input filename
    strict_map_key : bool, default=False
        If True, only allow string keys in msgpack dictionaries
        If False, allow integer and other keys in msgpack dictionaries
        
    Returns:
    --------
    atoms_list : list of Atoms
        The loaded atoms
    metadata : dict
        The metadata stored with the trajectory
    """
    if filename.endswith('.atraj'):
        return read_atraj(filename)
    elif filename.endswith('.tess'):
        return read_tess(filename)
    else:
        extn = filename.split('.')[-1]
        raise ValueError(f"Unsupported file extension: {extn}")
