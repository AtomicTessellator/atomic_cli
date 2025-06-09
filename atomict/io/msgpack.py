from typing import Union, List, Dict, Any, Tuple


def atoms_to_dict(atoms_list, selective=False):
    """Extract all properties from ASE Atoms objects into a standardized dictionary.
    
    Parameters:
    -----------
    atoms_list : List[Atoms]
        List of ASE Atoms objects
    selective : bool
        If True, only include non-default properties
        
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
    unique_symbols = set()
    for a in atoms_list:
        unique_symbols.update(a.get_chemical_symbols())
    unique_symbols = sorted(list(unique_symbols))
    
    # Store symbols data differently for variable atom count trajectories
    data['unique_symbols'] = unique_symbols
    data['symbols'] = []
    for a in atoms_list:
        # Convert each atom's symbols to indices in the unique_symbols list
        symbols_idx = [unique_symbols.index(s) for s in a.get_chemical_symbols()]
        data['symbols'].append(np.array(symbols_idx, dtype=np.uint16))
    
    # Store standard properties
    data['positions'] = [a.get_positions() for a in atoms_list]
    
    # Handle cell objects consistently
    cells = []
    for a in atoms_list:
        cell = a.get_cell()
        # Handle Cell object vs numpy array
        if hasattr(cell, 'array'):
            cells.append(np.array(cell.array, dtype=np.float32))
        else:
            cells.append(np.array(cell, dtype=np.float32))
    data['cell'] = cells
    
    data['pbc'] = [a.get_pbc() for a in atoms_list]
    data['numbers'] = [a.get_atomic_numbers() for a in atoms_list]
    
    # Always include masses for proper atomic weights
    data['masses'] = [a.get_masses() for a in atoms_list]
    
    # For selective mode, only include non-default properties
    if selective:
        # Include tags only if they're non-zero
        has_tags = any(np.any(a.get_tags() != 0) for a in atoms_list)
        if has_tags:
            data['tags'] = [a.get_tags() for a in atoms_list]
        
        # Include momenta only if they're non-zero
        has_momenta = any(np.any(np.abs(a.get_momenta()) > 1e-10) for a in atoms_list)
        if has_momenta:
            data['momenta'] = [a.get_momenta() for a in atoms_list]
        
        # Include charges only if they're non-zero
        has_charges = any(np.any(np.abs(a.get_initial_charges()) > 1e-10) for a in atoms_list)
        if has_charges:
            data['initial_charges'] = [a.get_initial_charges() for a in atoms_list]
        
        # Include magmoms only if they're non-zero
        has_magmoms = any(np.any(np.abs(a.get_initial_magnetic_moments()) > 1e-10) for a in atoms_list)
        if has_magmoms:
            data['initial_magmoms'] = [a.get_initial_magnetic_moments() for a in atoms_list]
    else:
        # Always include these for maximum compatibility
        data['tags'] = [a.get_tags() for a in atoms_list]
        data['momenta'] = [a.get_momenta() for a in atoms_list]
        data['initial_charges'] = [a.get_initial_charges() for a in atoms_list]
        data['initial_magmoms'] = [a.get_initial_magnetic_moments() for a in atoms_list]
    
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
        data['forces'] = [a.arrays.get('forces', np.zeros((len(a), 3), dtype=np.float32)) 
                           for a in atoms_list]
    
    # Handle calculator data - store in all cases where a calculator exists
    has_calc = False
    calc_data_list = []
    
    for a in atoms_list:
        calc_data = {}
        calc_found = False
        
        # First try getting data from the calculator object directly
        if hasattr(a, 'calc') and a.calc is not None:
            has_calc = True
            calc_found = True
            # Store calculator name and results
            calc_name = a.calc.__class__.__name__
            calc_data['name'] = calc_name
            
            # Try all standard properties
            for prop in ['energy', 'free_energy', 'forces', 'stress', 'dipole', 'charges', 'magmom', 'magmoms']:
                try:
                    if hasattr(a.calc, 'results') and prop in a.calc.results:
                        calc_data[prop] = a.calc.results[prop]
                    else:
                        value = a.calc.get_property(prop, a)
                        if value is not None:
                            calc_data[prop] = value
                except Exception:
                    pass
        
        # If no calculator directly available, try to get data from atoms.info
        if not calc_found and hasattr(a, 'info'):
            # Check for calculator data stored in info
            calc_name = a.info.get('_calc_name')
            
            if calc_name:
                has_calc = True
                calc_data['name'] = calc_name
                
                # Extract stored calculator properties
                for key, value in a.info.items():
                    if key.startswith('_calc_') and key != '_calc_name':
                        prop_name = key[6:]  # Remove '_calc_' prefix
                        calc_data[prop_name] = value
                
                # If we found any calculator info, mark as found
                if len(calc_data) > 1:  # More than just the name
                    calc_found = True
        
        calc_data_list.append(calc_data)
    
    if has_calc:
        data['calc_results'] = calc_data_list
    
    # Include stress only if present in any frame
    has_stress = any(hasattr(a, 'stress') and a.stress is not None for a in atoms_list)
    if has_stress:
        data['stress'] = [getattr(a, 'stress', np.zeros(6, dtype=np.float32)) for a in atoms_list]
    
    # Store atom info dictionaries
    if any(a.info for a in atoms_list):
        infos = []
        for a in atoms_list:
            info = a.info.copy()
            # Call to_dict on each info dictionary
            for key, value in info.items():
                if hasattr(value, 'to_dict') and callable(value.to_dict):
                    info[key] = value.to_dict()
                elif hasattr(value, 'todict') and callable(value.todict):
                    info[key] = value.todict()
                else:
                    info[key] = value
            infos.append(info)
        data['atom_infos'] = infos

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
        from ase.calculators.singlepoint import SinglePointCalculator
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")
    
    n_frames = data['n_frames']
    atoms_list = []
    
    # Get unique symbols
    unique_symbols = data['unique_symbols']
    symbols_map = data['symbols']
    
    # Loop through frames
    for i in range(n_frames):
        # Get symbols for this frame - handle both old and new format
        if isinstance(symbols_map[i], np.ndarray):
            frame_symbols = [unique_symbols[idx] for idx in symbols_map[i]]
        else:
            # Legacy format - symbols were stored as a 2D array
            idx = i * data['n_atoms'][i]
            frame_symbols = [unique_symbols[symbols_map[idx + j]] for j in range(data['n_atoms'][i])]
        
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
        
        if 'top_mask' in data and i < len(data['top_mask']) and data['top_mask'][i] is not None:
            atoms.top_mask = np.array(data['top_mask'][i], dtype=bool)

        if 'numbers' in data:
            atoms.set_atomic_numbers(data['numbers'][i])

        if 'constraints' in data and i < len(data['constraints']):
            for c in data['constraints'][i]:
                atoms.constraints.append(dict2constraint(c))

        if 'ase_objtype' in data and i < len(data['ase_objtype']) and data['ase_objtype'][i] is not None:
            atoms.ase_objtype = data['ase_objtype'][i]

        if 'forces' in data and i < len(data['forces']):
            atoms.arrays['forces'] = data['forces'][i]

        if 'stress' in data and i < len(data['stress']):
            atoms.stress = np.array(data['stress'][i], dtype=np.float64).copy()
        
        # Restore atom info
        if 'atom_infos' in data and i < len(data['atom_infos']):
            atoms.info.update(data['atom_infos'][i])
        
        # Restore custom arrays
        if 'custom_arrays' in data:
            for key, values in data['custom_arrays'].items():
                if i < len(values) and values[i] is not None:
                    atoms.arrays[key] = values[i]
        
        # Restore calculator if present
        calc_created = False
        calc_data = {}
        
        # First try from calc_results (new format)
        if 'calc_results' in data and i < len(data['calc_results']):
            calc_data = data['calc_results'][i]
            
            if calc_data and len(calc_data) > 1:  # Only create calculator if there's data beyond just the name
                # Initialize a SinglePointCalculator
                calc = SinglePointCalculator(atoms)
                
                # Set all available results directly to results dict
                for key, value in calc_data.items():
                    if key != 'name':  # Skip calculator name
                        calc.results[key] = value
                
                # Only set calculator if we have actual results
                if calc.results:
                    atoms.calc = calc
                    calc_created = True
        
        # If no calculator created yet, check atoms.info for calculator data
        if not calc_created:
            calc_info = {}
            for key, value in atoms.info.items():
                if key.startswith('_calc_') and key != '_calc_name':
                    prop_name = key[6:]  # Remove '_calc_' prefix
                    calc_info[prop_name] = value
            
            # Create calculator if we have any info data
            if calc_info:
                calc = SinglePointCalculator(atoms)
                for key, value in calc_info.items():
                    calc.results[key] = value
                atoms.calc = calc
        
        atoms_list.append(atoms)
    
    return atoms_list


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
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

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
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

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
    
    # Extract properties to dictionary - no selective mode for trajectories
    atoms_data = atoms_to_dict(atoms_list, selective=False)
    
    # Add atoms data to the trajectory container
    traj_data['atoms_data'] = atoms_data
    
    # Pack and save
    with open(filename, 'wb') as f:
        msgpack.pack(traj_data, f, use_bin_type=True)


def load_msgpack_trajectory(filename: str, strict_map_key: bool = True) -> Tuple[List['ase.Atoms'], Dict]:
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
    try:
        import msgpack
        import msgpack_numpy as m
    except ImportError:
        raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")

    # Enable numpy array deserialization
    m.patch()
    
    # Load data
    with open(filename, 'rb') as f:
        traj_data = msgpack.unpack(f, raw=False, strict_map_key=strict_map_key)
    
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