from typing import Union, List, Dict, Any, Tuple


def load_msgpack(filename: str) -> Union['ase.Atoms', List['ase.Atoms']]:
    """Load atoms from a msgpack file with high efficiency and speed."""

    try:
        import numpy as np
        import msgpack
        import msgpack_numpy as m
        from ase import Atoms
        from ase.constraints import dict2constraint
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

        # Always set the stress attribute, even if it's zeros
        if 'stress' in data:
            atoms.stress = np.array(data['stress'][i], dtype=np.float64).copy()
        else:
            # Set default stress if not in data
            atoms.stress = np.zeros(6, dtype=np.float64)

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
        # Check if cell is already a numpy array or if it's a Cell object
        if hasattr(cell, 'array'):
            cells.append(np.array(cell.array, dtype=np.float32))
        else:
            # Already a numpy array
            cells.append(np.array(cell, dtype=np.float32))
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
    
    # Include ase_objtype if it exists
    if any(hasattr(a, 'ase_objtype') for a in atoms_list):
        data['ase_objtype'] = [a.ase_objtype if hasattr(a, 'ase_objtype') else None for a in atoms_list]
    
    # Only include top_mask if it exists
    if any(hasattr(a, 'top_mask') for a in atoms_list):
        top_masks = []
        for a in atoms_list:
            if hasattr(a, 'top_mask'):
                top_masks.append(a.top_mask.tolist())
            else:
                top_masks.append(None)
        data['top_mask'] = top_masks
    
    # Include atomic numbers
    data['numbers'] = np.asarray([a.get_atomic_numbers() for a in atoms_list], dtype=np.int32)
    
    # Include constraints if they exist
    if any(a.constraints for a in atoms_list):
        data['constraints'] = [[c.todict() for c in a.constraints] for a in atoms_list]
    
    # Include forces if they exist
    if any('forces' in a.arrays for a in atoms_list):
        data['forces'] = np.asarray([a.arrays.get('forces', None) for a in atoms_list], dtype=np.float32)
    
    # Always include stress (zeros if not present)
    stresses = []
    for a in atoms_list:
        if hasattr(a, 'stress') and a.stress is not None:
            stresses.append(a.stress)
        else:
            stresses.append(np.zeros(6, dtype=np.float32))
    data['stress'] = np.asarray(stresses, dtype=np.float32)
    
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
    
    # Process each frame individually to support variable numbers of atoms per frame
    all_symbols = []
    symbols_per_frame = []
    positions_per_frame = []
    cells_per_frame = []
    pbc_per_frame = []
    tags_per_frame = []
    masses_per_frame = []
    momenta_per_frame = []
    charges_per_frame = []
    magmoms_per_frame = []
    
    # Process each frame to collect data
    for atoms in atoms_list:
        # Collect symbols
        symbols = atoms.get_chemical_symbols()
        all_symbols.extend(symbols)
        symbols_per_frame.append(symbols)
        
        # Collect other properties
        positions_per_frame.append(atoms.get_positions())
        
        # Handle the cell the same way as in save_msgpack
        cell = atoms.get_cell()
        if hasattr(cell, 'array'):
            cells_per_frame.append(np.array(cell.array, dtype=np.float32))
        else:
            # Already a numpy array
            cells_per_frame.append(np.array(cell, dtype=np.float32))
            
        pbc_per_frame.append(atoms.get_pbc())
        tags_per_frame.append(atoms.get_tags())
        masses_per_frame.append(atoms.get_masses())
        momenta_per_frame.append(atoms.get_momenta())
        charges_per_frame.append(atoms.get_initial_charges())
        magmoms_per_frame.append(atoms.get_initial_magnetic_moments())
    
    # Get unique symbols across all frames
    unique_symbols, inverse_map = np.unique(all_symbols, return_inverse=True)
    atoms_data['unique_symbols'] = unique_symbols.tolist()
    
    # Map the symbols for each frame separately
    start_idx = 0
    symbols_mapped = []
    for i, symbols in enumerate(symbols_per_frame):
        n_atoms = len(symbols)
        frame_map = inverse_map[start_idx:start_idx + n_atoms]
        symbols_mapped.append(frame_map.astype(np.uint16))
        start_idx += n_atoms
    
    atoms_data['symbols'] = symbols_mapped
    atoms_data['positions'] = [pos.astype(np.float32) for pos in positions_per_frame]
    atoms_data['cell'] = cells_per_frame
    atoms_data['pbc'] = pbc_per_frame
    
    # Always save masses to ensure custom masses are preserved
    atoms_data['masses'] = [masses.astype(np.float32) for masses in masses_per_frame]
    
    # Check if tags are non-default
    if any(np.any(tags != 0) for tags in tags_per_frame):
        atoms_data['tags'] = [tags.astype(np.int32) for tags in tags_per_frame]
    
    # Check if momenta are non-zero
    if any(np.any(mom) for mom in momenta_per_frame):
        atoms_data['momenta'] = [mom.astype(np.float32) for mom in momenta_per_frame]
    
    # Check if initial charges are non-zero
    if any(np.any(chg) for chg in charges_per_frame):
        atoms_data['initial_charges'] = [chg.astype(np.float32) for chg in charges_per_frame]
    
    # Check if magnetic moments are non-zero
    if any(np.any(mag) for mag in magmoms_per_frame):
        atoms_data['initial_magmoms'] = [mag.astype(np.float32) for mag in magmoms_per_frame]
    
    # Store ALL atom-specific info dictionaries
    atom_infos = []
    for atom in atoms_list:
        # Store complete info dictionary
        atom_infos.append(atom.info.copy())
    
    if any(atom_infos):
        atoms_data['atom_infos'] = atom_infos
    
    # Store custom arrays (any arrays not in the standard set)
    standard_arrays = {'numbers', 'positions', 'momenta', 'masses', 'tags', 'charges'}
    custom_arrays_by_frame = []
    
    for atom in atoms_list:
        custom_arrays = {}
        for key, value in atom.arrays.items():
            if key not in standard_arrays:
                custom_arrays[key] = value
        custom_arrays_by_frame.append(custom_arrays)
    
    # Only store if there are custom arrays
    if any(custom_arrays_by_frame):
        atoms_data['custom_arrays'] = custom_arrays_by_frame
    
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
    
    # Handle both the old version (array for all frames) and new version (list of arrays per frame)
    symbols_map = atoms_data['symbols']
    
    # Loop through frames
    for i in range(n_frames):
        # Check if symbols_map is a list of arrays (new version) or a single array (old version)
        if isinstance(symbols_map, list):
            # New version - each frame has its own symbols map
            frame_symbols = [unique_symbols[idx] for idx in symbols_map[i]]
        else:
            # Old version - reshape the global map
            atoms_per_frame = atoms_data['n_atoms'][i]
            start_idx = sum(atoms_data['n_atoms'][:i])
            end_idx = start_idx + atoms_per_frame
            frame_symbols = [unique_symbols[idx] for idx in symbols_map[start_idx:end_idx]]
        
        # Positions and cell may be a list of arrays (new version) or a single array (old version)
        positions = atoms_data['positions'][i] if isinstance(atoms_data['positions'], list) else atoms_data['positions'][i]
        cell = atoms_data['cell'][i] if isinstance(atoms_data['cell'], list) else atoms_data['cell'][i]
        pbc = atoms_data['pbc'][i] if isinstance(atoms_data['pbc'], list) else atoms_data['pbc'][i]
        
        # Create atoms object
        atoms = Atoms(
            symbols=frame_symbols,
            positions=positions,
            cell=cell,
            pbc=pbc,
        )
        
        # Set optional properties if they exist
        # Handle both list-of-arrays and single-array versions
        if 'tags' in atoms_data:
            tags = atoms_data['tags'][i] if isinstance(atoms_data['tags'], list) else atoms_data['tags'][i]
            atoms.set_tags(tags)
        
        if 'masses' in atoms_data:
            masses = atoms_data['masses'][i] if isinstance(atoms_data['masses'], list) else atoms_data['masses'][i]
            atoms.set_masses(masses)
        
        if 'momenta' in atoms_data:
            momenta = atoms_data['momenta'][i] if isinstance(atoms_data['momenta'], list) else atoms_data['momenta'][i]
            atoms.set_momenta(momenta)
        
        if 'initial_charges' in atoms_data:
            charges = atoms_data['initial_charges'][i] if isinstance(atoms_data['initial_charges'], list) else atoms_data['initial_charges'][i]
            atoms.set_initial_charges(charges)
        
        if 'initial_magmoms' in atoms_data:
            magmoms = atoms_data['initial_magmoms'][i] if isinstance(atoms_data['initial_magmoms'], list) else atoms_data['initial_magmoms'][i]
            atoms.set_initial_magnetic_moments(magmoms)
        
        # Restore atom-specific info
        if 'atom_infos' in atoms_data and i < len(atoms_data['atom_infos']):
            atoms.info.update(atoms_data['atom_infos'][i])
        
        # Restore custom arrays
        if 'custom_arrays' in atoms_data and i < len(atoms_data['custom_arrays']):
            custom_arrays = atoms_data['custom_arrays'][i]
            for key, value in custom_arrays.items():
                atoms.arrays[key] = value
        
        atoms_list.append(atoms)
    
    return atoms_list, metadata