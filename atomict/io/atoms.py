from typing import List
import numpy as np
from ase import Atoms
from ase.constraints import dict2constraint
from ase.calculators.singlepoint import SinglePointCalculator


def atoms_to_dict(atoms_list: List[Atoms], selective: bool = False):
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
    
    n_frames = data['n_frames']
    atoms_list = [None] * n_frames  # type: ignore[assignment]

    # Cache lookups outside loop for speed
    unique_symbols = data['unique_symbols']
    unique_symbols_array = np.asarray(unique_symbols, dtype=object)
    symbols_map = data['symbols']
    positions = data['positions']
    cells = data['cell']
    pbc = data['pbc']
    n_atoms_seq = data['n_atoms']

    numbers_data = data.get('numbers')
    has_numbers = numbers_data is not None

    tags = data.get('tags')
    masses = data.get('masses')
    momenta = data.get('momenta')
    initial_charges = data.get('initial_charges')
    initial_magmoms = data.get('initial_magmoms')
    top_mask = data.get('top_mask')
    constraints = data.get('constraints')
    ase_objtype = data.get('ase_objtype')
    forces = data.get('forces')
    stress = data.get('stress')
    atom_infos = data.get('atom_infos')
    custom_arrays = data.get('custom_arrays')
    calc_results = data.get('calc_results')

    for i in range(n_frames):
        symbols_entry = symbols_map[i]
        if isinstance(symbols_entry, np.ndarray):
            symbols_idx = symbols_entry
        else:
            start_idx = i * n_atoms_seq[i]
            symbols_idx = symbols_map[start_idx:start_idx + n_atoms_seq[i]]

        if has_numbers:
            atoms = Atoms(
                numbers=numbers_data[i],
                positions=positions[i],
                cell=cells[i],
                pbc=pbc[i],
            )
        else:
            # Avoid numpy conversion and tolist() for small arrays
            if isinstance(symbols_idx, np.ndarray):
                frame_symbols = [unique_symbols[idx] for idx in symbols_idx]
            else:
                frame_symbols = [unique_symbols[idx] for idx in symbols_idx]
            atoms = Atoms(
                symbols=frame_symbols,
                positions=positions[i],
                cell=cells[i],
                pbc=pbc[i],
            )
        
        # Set optional properties if they exist
        if tags is not None:
            atoms.set_tags(tags[i])

        if masses is not None:
            atoms.set_masses(masses[i])

        if momenta is not None:
            atoms.set_momenta(momenta[i])

        if initial_charges is not None:
            atoms.set_initial_charges(initial_charges[i])

        if initial_magmoms is not None:
            atoms.set_initial_magnetic_moments(initial_magmoms[i])

        if top_mask is not None and i < len(top_mask) and top_mask[i] is not None:
            atoms.top_mask = np.asarray(top_mask[i], dtype=bool)

        if constraints is not None and i < len(constraints):
            for c in constraints[i]:
                atoms.constraints.append(dict2constraint(c))

        if ase_objtype is not None and i < len(ase_objtype) and ase_objtype[i] is not None:
            atoms.ase_objtype = ase_objtype[i]

        if forces is not None and i < len(forces):
            atoms.arrays['forces'] = forces[i]

        if stress is not None and i < len(stress):
            atoms.stress = np.array(stress[i], dtype=np.float64).copy()
        
        # Restore atom info
        if atom_infos is not None and i < len(atom_infos):
            atoms.info.update(atom_infos[i])
        
        # Restore custom arrays
        if custom_arrays is not None:
            for key, values in custom_arrays.items():
                if i < len(values) and values[i] is not None:
                    atoms.arrays[key] = values[i]
        
        # Restore calculator if present
        calc_created = False
        calc_data = {}

        # First try from calc_results (new format)
        if calc_results is not None and i < len(calc_results):
            calc_data = calc_results[i]
            
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
        
        atoms_list[i] = atoms
    
    return atoms_list
