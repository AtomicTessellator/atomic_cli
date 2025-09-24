from typing import Union, List, Dict, Tuple
from atomict.io.formats.atraj import write_atraj, read_atraj
from atomict.io.formats.tess import write_tess, read_tess
from atomict.io.atoms import dict_to_atoms, atoms_to_dict


def load_msgpack(filename: str, strict_map_key: bool = True) -> Union['ase.Atoms', List['ase.Atoms']]:
    """Load atoms from a msgpack file with high efficiency and speed.
    
    Parameters:
    -----------
    filename : str
        The input filename
    strict_map_key : bool, default=False
        If True, only allow string keys in msgpack dictionaries
        If False, allow integer and other keys in msgpack dictionaries
    """

    try:
        import msgpack
        import msgpack_numpy as m
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()
    
    # Load data
    with open(filename, 'rb') as f:
        data = msgpack.unpack(f, raw=False, strict_map_key=strict_map_key)
    
    # Convert to atoms objects
    atoms_list = dict_to_atoms(data)
    
    # Return single atom or list based on input
    return atoms_list[0] if data['n_frames'] == 1 else atoms_list


def save_msgpack(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str):
    """Save atoms to a msgpack file with high efficiency and speed."""

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[utils]` to use msgpack I/O")

    # Enable numpy array serialization
    m.patch()
    
    # Single atoms case - convert to list
    if isinstance(atoms, Atoms):
        atoms_list = [atoms]
    else:
        atoms_list = atoms
    
    # Extract properties to dictionary - use selective mode for single atoms
    # to avoid storing default properties
    selective = len(atoms_list) == 1
    data = atoms_to_dict(atoms_list, selective=selective)
    
    # Pack and save
    with open(filename, 'wb') as f:
        msgpack.pack(data, f, use_bin_type=True)


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
