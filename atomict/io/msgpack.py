from typing import Union, List


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
    data['cell'] = np.asarray([np.array(a.get_cell()) for a in atoms_list], dtype=np.float32)
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