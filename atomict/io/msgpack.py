from typing import Union, List, Dict, Any, Tuple


def atoms_to_dict(atoms_list):
    """Extract all properties from ASE Atoms objects into a standardized dictionary.
    
    Parameters:
    -----------
    atoms_list : List[Atoms]
        List of ASE Atoms objects
        
    Returns:
    --------
    Dict
        Dictionary with all extracted properties
    """
    import numpy as np
    
    # Create data structure with common properties
    data = {
        'n_frames': len(atoms_list),
        'n_atoms': [len(a) for a in atoms_list],
    }
    
    # Process all symbols efficiently
    all_symbols = []
    for a in atoms_list:
        all_symbols.extend(a.get_chemical_symbols())
    unique_symbols, symbols_map = np.unique(all_symbols, return_inverse=True)
    
    # Store symbols data
    data['unique_symbols'] = unique_symbols.tolist()
    data['symbols'] = symbols_map.reshape([len(atoms_list), -1]).astype(np.uint16)
    
    # Store standard properties
    data['positions'] = np.asarray([a.get_positions() for a in atoms_list], dtype=np.float32)
    
    # Handle cell objects consistently
    cells = []
    for a in atoms_list:
        cell = a.get_cell()
        # Handle Cell object vs numpy array
        if hasattr(cell, 'array'):
            cells.append(np.array(cell.array, dtype=np.float32))
        else:
            cells.append(np.array(cell, dtype=np.float32))
    data['cell'] = np.array(cells, dtype=np.float32)
    
    data['pbc'] = np.asarray([a.get_pbc() for a in atoms_list], dtype=bool)
    data['numbers'] = np.asarray([a.get_atomic_numbers() for a in atoms_list], dtype=np.int32)
    
    # Always include these properties for consistency
    data['tags'] = np.asarray([a.get_tags() for a in atoms_list], dtype=np.int32)
    data['masses'] = np.asarray([a.get_masses() for a in atoms_list], dtype=np.float32)
    data['momenta'] = np.asarray([a.get_momenta() for a in atoms_list], dtype=np.float32)
    data['initial_charges'] = np.asarray([a.get_initial_charges() for a in atoms_list], dtype=np.float32)
    data['initial_magmoms'] = np.asarray([a.get_initial_magnetic_moments() for a in atoms_list], dtype=np.float32)
    
    # Get all constraints
    if any(a.constraints for a in atoms_list):
        data['constraints'] = [[c.todict() for c in a.constraints] for a in atoms_list]
    
    # Handle custom properties
    if any(hasattr(a, 'ase_objtype') for a in atoms_list):
        data['ase_objtype'] = [getattr(a, 'ase_objtype', None) for a in atoms_list]
    
    if any(hasattr(a, 'top_mask') for a in atoms_list):
        data['top_mask'] = [getattr(a, 'top_mask', None) for a in atoms_list]
    
    # Handle forces array
    if any('forces' in a.arrays for a in atoms_list):
        data['forces'] = np.asarray([a.arrays.get('forces', np.zeros((len(a), 3), dtype=np.float32)) 
                                     for a in atoms_list], dtype=np.float32)
    
    # Always include stress (zeros if not present)
    data['stress'] = np.asarray([getattr(a, 'stress', np.zeros(6, dtype=np.float32)) 
                                for a in atoms_list], dtype=np.float32)
    
    # Store atom info dictionaries
    if any(a.info for a in atoms_list):
        data['atom_infos'] = [a.info.copy() for a in atoms_list]
    
    # Extract custom arrays
    standard_arrays = {'numbers', 'positions', 'momenta', 'masses', 'tags', 'charges'}
    custom_arrays = {}
    
    for i, atom in enumerate(atoms_list):
        for key, value in atom.arrays.items():
            if key not in standard_arrays:
                if key not in custom_arrays:
                    custom_arrays[key] = [None] * len(atoms_list)
                custom_arrays[key][i] = value
    
    if custom_arrays:
        data['custom_arrays'] = custom_arrays
    
    return data


def dict_to_atoms(data):
    """Create ASE Atoms objects from a dictionary of properties.
    
    Parameters:
    -----------
    data : Dict
        Dictionary with all properties
        
    Returns:
    --------
    List[Atoms]
        List of ASE Atoms objects
    """
    try:
        import numpy as np
        from ase import Atoms
        from ase.constraints import dict2constraint
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")
    
    n_frames = data['n_frames']
    atoms_list = []
    
    # Get unique symbols
    unique_symbols = data['unique_symbols']
    symbols_map = data['symbols']
    
    # Loop through frames
    for i in range(n_frames):
        # Get symbols for this frame
        frame_symbols = [unique_symbols[idx] for idx in symbols_map[i]]
        
        # Create atoms object with basic properties
        atoms = Atoms(
            symbols=frame_symbols,
            positions=data['positions'][i],
            cell=data['cell'][i],
            pbc=data['pbc'][i],
        )
        
        # Set optional properties if they exist
        if 'tags' in data:
            atoms.set_tags(data['tags'][i])
        
        if 'masses' in data:
            atoms.set_masses(data['masses'][i])
        
        if 'momenta' in data:
            atoms.set_momenta(data['momenta'][i])
        
        if 'initial_charges' in data:
            atoms.set_initial_charges(data['initial_charges'][i])
        
        if 'initial_magmoms' in data:
            atoms.set_initial_magnetic_moments(data['initial_magmoms'][i])
        
        if 'top_mask' in data and data['top_mask'][i] is not None:
            atoms.top_mask = np.array(data['top_mask'][i], dtype=bool)

        if 'numbers' in data:
            atoms.set_atomic_numbers(data['numbers'][i])

        if 'constraints' in data:
            for c in data['constraints'][i]:
                atoms.constraints.append(dict2constraint(c))

        if 'ase_objtype' in data and data['ase_objtype'][i] is not None:
            atoms.ase_objtype = data['ase_objtype'][i]

        if 'forces' in data:
            atoms.arrays['forces'] = data['forces'][i]

        if 'stress' in data:
            atoms.stress = np.array(data['stress'][i], dtype=np.float64).copy()
        
        # Restore atom info
        if 'atom_infos' in data and i < len(data['atom_infos']):
            atoms.info.update(data['atom_infos'][i])
        
        # Restore custom arrays
        if 'custom_arrays' in data:
            for key, values in data['custom_arrays'].items():
                if values[i] is not None:
                    atoms.arrays[key] = values[i]
        
        atoms_list.append(atoms)
    
    return atoms_list


def load_msgpack(filename: str) -> Union['ase.Atoms', List['ase.Atoms']]:
    """Load atoms from a msgpack file with high efficiency and speed."""

    try:
        import msgpack
        import msgpack_numpy as m
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()
    
    # Load data
    with open(filename, 'rb') as f:
        data = msgpack.unpack(f, raw=False)
    
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
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

    # Enable numpy array serialization
    m.patch()
    
    # Single atoms case - convert to list
    if isinstance(atoms, Atoms):
        atoms_list = [atoms]
    else:
        atoms_list = atoms
    
    # Extract properties to dictionary
    data = atoms_to_dict(atoms_list)
    
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
    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

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
    
    # Extract properties to dictionary
    atoms_data = atoms_to_dict(atoms_list)
    
    # Add atoms data to the trajectory container
    traj_data['atoms_data'] = atoms_data
    
    # Pack and save
    with open(filename, 'wb') as f:
        msgpack.pack(traj_data, f, use_bin_type=True)


def load_msgpack_trajectory(filename: str) -> Tuple[List['ase.Atoms'], Dict]:
    """Load atoms from a msgpack trajectory file with metadata.
    
    Parameters:
    -----------
    filename : str
        The input filename
        
    Returns:
    --------
    atoms_list : list of Atoms
        The loaded atoms
    metadata : dict
        The metadata stored with the trajectory
    """
    try:
        import msgpack
        import msgpack_numpy as m
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()
    
    # Load data
    with open(filename, 'rb') as f:
        traj_data = msgpack.unpack(f, raw=False)
    
    # Check if this is a new-style trajectory with format_version
    if isinstance(traj_data, dict) and 'format_version' in traj_data:
        metadata = traj_data.get('metadata', {})
        atoms_data = traj_data.get('atoms_data', {})
    else:
        # Legacy format - just raw atoms data
        metadata = {}
        atoms_data = traj_data
    
    # Convert to atoms objects
    atoms_list = dict_to_atoms(atoms_data)
    
    # Make sure atoms_list is always a list
    if not isinstance(atoms_list, list):
        atoms_list = [atoms_list]
    
    return atoms_list, metadata