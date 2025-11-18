from typing import Union, List, Dict, TYPE_CHECKING


if TYPE_CHECKING:
    import ase


def write_atraj(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str, metadata: Dict = None):
    """Save atoms to a mpv1 file."""

    from atomict.io.atoms import atoms_to_dict
    import msgpack
    import msgpack_numpy as m
    from ase import Atoms

    # Enable numpy array serialization
    m.patch()
    
    # Single atoms case - convert to list
    if isinstance(atoms, Atoms):
        atoms_list = [atoms]
    else:
        atoms_list = atoms
    
    # Create container for the trajectory data
    traj_data = {
        'format_version': 1,  # Version for future compatibility
        'metadata': metadata or {},
    }
    
    # Extract properties to dictionary - no selective mode for trajectories
    atoms_data = atoms_to_dict(atoms_list, selective=False)
    
    # Add atoms data to the trajectory container
    traj_data['atoms_data'] = atoms_data
    
    # Pack and save
    with open(filename, 'wb') as f:
        msgpack.pack(traj_data, f, use_bin_type=True)


def read_atraj(filename: str) -> Union['ase.Atoms', List['ase.Atoms']]:

    """Load atoms from a mpv1 file."""
    from atomict.io.atoms import dict_to_atoms
    import msgpack
    import msgpack_numpy as m

    # Enable numpy array deserialization
    m.patch()
    
    # Load data
    with open(filename, 'rb') as f:
        traj_data = msgpack.unpack(f, raw=False, strict_map_key=True)
    
    # Check if this is a new-style trajectory with format_version
    if isinstance(traj_data, dict) and 'format_version' in traj_data:
        metadata = traj_data.get('metadata', {})
        atoms_data = traj_data.get('atoms_data', {})
    else:
        # Legacy format - just raw atoms data
        metadata = {}
        atoms_data = traj_data
    
    # Ensure that calculated properties are transferred to the calculator in dict_to_atoms
    if 'calc_results' not in atoms_data and hasattr(atoms_data, 'get') and atoms_data.get('forces') is not None:
        # If we have forces in the data but no calc_results, create calc_results entries
        calc_data_list = []
        n_frames = atoms_data.get('n_frames', 0)
        
        for i in range(n_frames):
            calc_data = {'name': 'SinglePointCalculator'}
            if 'forces' in atoms_data and i < len(atoms_data['forces']):
                calc_data['forces'] = atoms_data['forces'][i]
            if 'stress' in atoms_data and i < len(atoms_data['stress']):
                calc_data['stress'] = atoms_data['stress'][i]
            if 'energy' in atoms_data and i < len(atoms_data['energy']):
                calc_data['energy'] = atoms_data['energy'][i]
            calc_data_list.append(calc_data)
        
        atoms_data['calc_results'] = calc_data_list
    
    # Convert to atoms objects
    atoms_list = dict_to_atoms(atoms_data)
    
    # Make sure atoms_list is always a list
    if not isinstance(atoms_list, list):
        atoms_list = [atoms_list]
    
    return atoms_list, metadata

