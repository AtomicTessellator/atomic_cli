from typing import Union, List, Dict, Any, Tuple


def load_msgpack(filename: str) -> Union['ase.Atoms', List['ase.Atoms']]:
    """Load atoms from a msgpack file with high efficiency and speed."""

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()
    
    # Load data
    with open(filename, 'rb') as f:
        data = msgpack.unpack(f, raw=False)
    
    n_frames = data['n_frames']
    atoms_list = []
    
    # Get unique symbols
    unique_symbols = data['unique_symbols']
    symbols_map = data['symbols']
    
    # Loop through frames
    for i in range(n_frames):
        # Get symbols for this frame
        frame_symbols = [unique_symbols[idx] for idx in symbols_map[i]]
        
        # Create atoms object
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
        
        atoms_list.append(atoms)
    
    # Return single atom or list based on input
    return atoms_list[0] if n_frames == 1 else atoms_list


def save_msgpack(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str):
    """Save atoms to a msgpack file with high efficiency and speed."""

    try:
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
        import numpy as np
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

    # Enable numpy array serialization
    m.patch()
    
    # Single atoms case - convert to list
    if isinstance(atoms, Atoms):
        atoms_list = [atoms]
    else:
        atoms_list = atoms
    
    # Create data structure optimized for msgpack
    data = {
        'n_frames': len(atoms_list),
        'n_atoms': [len(a) for a in atoms_list],
    }
    
    # Collect data optimally
    # Store symbols as integers for efficiency
    all_symbols = []
    for a in atoms_list:
        all_symbols.extend(a.get_chemical_symbols())
    unique_symbols, symbols_map = np.unique(all_symbols, return_inverse=True)
    
    # Store data efficiently
    data['unique_symbols'] = unique_symbols.tolist()
    data['symbols'] = symbols_map.reshape([len(atoms_list), -1]).astype(np.uint16)
    
    # Store positions as float32 for better space efficiency
    data['positions'] = np.asarray([a.get_positions() for a in atoms_list], dtype=np.float32)
    
    # Fix for NumPy 2.0+ compatibility - explicitly convert Cell to array
    cells = []
    for a in atoms_list:
        cell = a.get_cell()
        cells.append(np.array(cell.array, dtype=np.float32))
    data['cell'] = np.array(cells, dtype=np.float32)
    
    data['pbc'] = np.asarray([a.get_pbc() for a in atoms_list], dtype=bool)
    
    # Only include non-default properties if they have values
    # Check first atom to see if we need to include these properties
    if any(atoms_list[0].get_tags() != 0):
        data['tags'] = np.asarray([a.get_tags() for a in atoms_list], dtype=np.int32)
    
    # Check if masses are non-default
    default_masses = atoms_list[0].get_masses() / atoms_list[0].get_atomic_numbers()
    if not np.allclose(default_masses, default_masses[0], rtol=1e-5):
        data['masses'] = np.asarray([a.get_masses() for a in atoms_list], dtype=np.float32)
    
    # Only include momenta if non-zero
    if np.any([np.any(a.get_momenta()) for a in atoms_list]):
        data['momenta'] = np.asarray([a.get_momenta() for a in atoms_list], dtype=np.float32)
    
    # Only include charges if non-zero
    if np.any([np.any(a.get_initial_charges()) for a in atoms_list]):
        data['initial_charges'] = np.asarray([a.get_initial_charges() for a in atoms_list], dtype=np.float32)
    
    # Only include magnetic moments if non-zero
    if np.any([np.any(a.get_initial_magnetic_moments()) for a in atoms_list]):
        data['initial_magmoms'] = np.asarray([a.get_initial_magnetic_moments() for a in atoms_list], dtype=np.float32)
    
    # Pack and save
    with open(filename, 'wb') as f:
        msgpack.pack(data, f, use_bin_type=True)


def save_msgpack_trajectory(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str, metadata: Dict = None):
    """Save atoms to a msgpack trajectory file with metadata.
    
    This version includes special handling for trajectory metadata like descriptions,
    calculator data, and constraints.
    
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
        import numpy as np
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
        'n_frames': len(atoms_list),
    }
    
    # Create data structure optimized for msgpack for the atoms data
    atoms_data = {
        'n_frames': len(atoms_list),
        'n_atoms': [len(a) for a in atoms_list],
    }
    
    # Collect data optimally
    # Store symbols as integers for efficiency
    all_symbols = []
    for a in atoms_list:
        all_symbols.extend(a.get_chemical_symbols())
    unique_symbols, symbols_map = np.unique(all_symbols, return_inverse=True)
    
    # Store data efficiently
    atoms_data['unique_symbols'] = unique_symbols.tolist()
    atoms_data['symbols'] = symbols_map.reshape([len(atoms_list), -1]).astype(np.uint16)
    
    # Store positions as float32 for better space efficiency
    atoms_data['positions'] = np.asarray([a.get_positions() for a in atoms_list], dtype=np.float32)
    
    # Fix for NumPy 2.0+ compatibility - explicitly convert Cell to array
    cells = []
    for a in atoms_list:
        cell = a.get_cell()
        cells.append(np.array(cell.array, dtype=np.float32))
    atoms_data['cell'] = np.array(cells, dtype=np.float32)
    
    atoms_data['pbc'] = np.asarray([a.get_pbc() for a in atoms_list], dtype=bool)
    
    # Only include non-default properties if they have values
    # Check first atom to see if we need to include these properties
    if any(atoms_list[0].get_tags() != 0):
        atoms_data['tags'] = np.asarray([a.get_tags() for a in atoms_list], dtype=np.int32)
    
    # Check if masses are non-default
    default_masses = atoms_list[0].get_masses() / atoms_list[0].get_atomic_numbers()
    if not np.allclose(default_masses, default_masses[0], rtol=1e-5):
        atoms_data['masses'] = np.asarray([a.get_masses() for a in atoms_list], dtype=np.float32)
    
    # Only include momenta if non-zero
    if np.any([np.any(a.get_momenta()) for a in atoms_list]):
        atoms_data['momenta'] = np.asarray([a.get_momenta() for a in atoms_list], dtype=np.float32)
    
    # Only include charges if non-zero
    if np.any([np.any(a.get_initial_charges()) for a in atoms_list]):
        atoms_data['initial_charges'] = np.asarray([a.get_initial_charges() for a in atoms_list], dtype=np.float32)
    
    # Only include magnetic moments if non-zero
    if np.any([np.any(a.get_initial_magnetic_moments()) for a in atoms_list]):
        atoms_data['initial_magmoms'] = np.asarray([a.get_initial_magnetic_moments() for a in atoms_list], dtype=np.float32)
    
    # Store atom-specific info dictionaries (for calculator data, constraints, etc.)
    atom_infos = []
    for atom in atoms_list:
        info = {}
        for key, value in atom.info.items():
            # Special keys that need special handling for serialization
            if key in ['_calc_data', '_constraints', '_traj_description', '_description', '_ase_version']:
                info[key] = value
        atom_infos.append(info)
    
    if any(atom_infos):
        atoms_data['atom_infos'] = atom_infos
    
    # Add atoms data to the trajectory container
    traj_data['atoms_data'] = atoms_data
    
    # Pack and save
    with open(filename, 'wb') as f:
        msgpack.pack(traj_data, f, use_bin_type=True)


def load_msgpack_trajectory(filename: str) -> Tuple[List['ase.Atoms'], Dict]:
    """Load atoms from a msgpack trajectory file with metadata.
    
    This version includes special handling for trajectory metadata.
    
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
        from ase import Atoms
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
    
    n_frames = atoms_data['n_frames']
    atoms_list = []
    
    # Get unique symbols
    unique_symbols = atoms_data['unique_symbols']
    symbols_map = atoms_data['symbols']
    
    # Loop through frames
    for i in range(n_frames):
        # Get symbols for this frame
        frame_symbols = [unique_symbols[idx] for idx in symbols_map[i]]
        
        # Create atoms object
        atoms = Atoms(
            symbols=frame_symbols,
            positions=atoms_data['positions'][i],
            cell=atoms_data['cell'][i],
            pbc=atoms_data['pbc'][i],
        )
        
        # Set optional properties if they exist
        if 'tags' in atoms_data:
            atoms.set_tags(atoms_data['tags'][i])
        
        if 'masses' in atoms_data:
            atoms.set_masses(atoms_data['masses'][i])
        
        if 'momenta' in atoms_data:
            atoms.set_momenta(atoms_data['momenta'][i])
        
        if 'initial_charges' in atoms_data:
            atoms.set_initial_charges(atoms_data['initial_charges'][i])
        
        if 'initial_magmoms' in atoms_data:
            atoms.set_initial_magnetic_moments(atoms_data['initial_magmoms'][i])
        
        # Restore atom-specific info
        if 'atom_infos' in atoms_data:
            atoms.info.update(atoms_data['atom_infos'][i])
        
        atoms_list.append(atoms)
    
    return atoms_list, metadata