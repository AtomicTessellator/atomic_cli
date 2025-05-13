import pytest
import numpy as np
from ase.build import molecule, bulk
from ase import Atom
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from atomict.ase.atatoms import ATAtoms


def test_basic_initialization():
    """Test basic initialization of ATAtoms."""
    atoms_orig = molecule('H2O')
    atoms = ATAtoms(atoms_orig)
    assert len(atoms) == 3


def test_atomic_position_manipulation():
    """Test manipulating atomic positions."""
    atoms = ATAtoms(molecule('H2O'))
    
    # Test direct position modification
    orig_pos = atoms[0].position.copy()
    atoms[0].position += [0.1, 0.2, 0.3]
    assert np.allclose(atoms[0].position, orig_pos + [0.1, 0.2, 0.3])
    
    # Test position indexing and assignment
    atoms[1].position[0] = 1.5
    assert np.isclose(atoms[1].position[0], 1.5)


def test_structural_transformations():
    """Test structural transformations like translation and rotation."""
    atoms = ATAtoms(molecule('CH4'))
    orig_pos = atoms.get_positions().copy()
    
    # Test translation
    atoms.translate([1.0, 0.0, 0.0])
    assert np.allclose(atoms.get_positions(), orig_pos + [1.0, 0.0, 0.0])
    
    # Test rotation
    atoms.rotate(90, 'z')
    # Check that rotation happened (positions changed)
    assert not np.allclose(atoms.get_positions(), orig_pos + [1.0, 0.0, 0.0])


def test_cell_operations():
    """Test operations on the cell."""
    atoms = ATAtoms(bulk('Cu'))
    orig_cell = atoms.get_cell().array.copy()
    
    # Test cell scaling
    atoms.set_cell(atoms.get_cell() * 1.05, scale_atoms=True)
    assert np.allclose(atoms.get_cell().array, orig_cell * 1.05)
    
    # Test adding vacuum
    atoms = ATAtoms(molecule('CO2'))
    pre_vacuum_volume = atoms.get_cell().volume
    atoms.center(vacuum=5.0)
    assert atoms.get_cell().volume > pre_vacuum_volume


def test_atom_addition_deletion():
    """Test adding and removing atoms."""
    atoms = ATAtoms(molecule('H2'))
    n_atoms = len(atoms)
    
    # Add an atom
    atoms += Atom('O', position=(0, 0, 0))
    assert len(atoms) == n_atoms + 1
    assert atoms.get_chemical_symbols()[-1] == 'O'
    
    # Delete an atom
    del atoms[0]
    assert len(atoms) == n_atoms


def test_slice_operations():
    """Test slicing operations on ATAtoms."""
    atoms = ATAtoms(molecule('CH4'))
    sliced = atoms[1:3]
    assert len(sliced) == 2
    assert isinstance(sliced, ATAtoms)


def test_repeat_operation():
    """Test the repeat operation."""
    atoms = ATAtoms(bulk('Al'))
    n_atoms = len(atoms)
    repeated = atoms.repeat((2, 2, 2))
    assert len(repeated) == n_atoms * 8
    assert isinstance(repeated, ATAtoms)


def test_property_setting():
    """Test setting various properties on atoms."""
    atoms = ATAtoms(molecule('O2'))
    
    # Test setting charges
    atoms.set_initial_charges([0.5, -0.5])
    assert np.isclose(atoms[0].charge, 0.5)
    assert np.isclose(atoms[1].charge, -0.5)
    
    # Test setting magnetic moments
    atoms.set_initial_magnetic_moments([1.0, -1.0])
    assert np.isclose(atoms.get_initial_magnetic_moments()[0], 1.0)
    assert np.isclose(atoms.get_initial_magnetic_moments()[1], -1.0)
    
    # Test custom arrays
    atoms.set_array('custom_data', np.array([3.14, 2.71]))
    assert np.isclose(atoms.get_array('custom_data')[0], 3.14)


def test_positions_assignment():
    """Test direct assignment of positions."""
    atoms = ATAtoms(molecule('H2'))
    new_positions = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    atoms.set_positions(new_positions)
    assert np.allclose(atoms.positions, new_positions)


# def test_delta_recording():
#     """Test recording of deltas during operations."""
#     atoms = ATAtoms(molecule('H2'), batch_size=100)  # Avoid actual syncing
#     initial_diff_ct = len(atoms._diffs)
    
#     # Perform several operations
#     atoms.translate([1.0, 0.0, 0.0])
#     atoms.rotate(45, 'z')
#     atoms[0].position += [0.5, 0.0, 0.0]
    
#     assert len(atoms._diffs) > initial_diff_ct


# def test_optimization_workflow():
#     """Test optimization workflow with a calculator."""
#     atoms = ATAtoms(bulk('Cu', cubic=True))
#     atoms.rattle(stdev=0.1)  # Displace atoms randomly to create forces
    
#     # Save initial positions and delta count
#     initial_positions = atoms.positions.copy()
#     initial_delta_count = len(atoms._diffs)
    
#     # Setup calculator and optimizer
#     atoms.calc = EMT()
#     opt = BFGS(atoms)
    
#     # Run a few steps of optimization
#     opt.run(fmax=0.01, steps=3)
    
#     # Check that positions changed
#     assert not np.allclose(atoms.positions, initial_positions)
    
#     # Check that we recorded deltas during optimization
#     assert len(atoms._diffs) > initial_delta_count
    
#     # Check that forces were lowered
#     forces = atoms.get_forces()
#     max_force = np.max(np.abs(forces))
#     assert max_force < 1.0  # Using a relaxed criterion


# def test_deepdiff_state_tracking():
#     """Test DeepDiff state tracking."""
#     # Create a new system and make a change to ensure we have diffs
#     atoms = ATAtoms(molecule('H2O'), batch_size=100)
    
#     # Verify initial state is recorded
#     assert len(atoms._diffs) >= 1
#     assert 'state' in atoms._diffs[0]
    
#     # Create state changes by translating
#     atoms.translate([1.0, 0.0, 0.0])
    
#     # Check if we have more than one diff after making changes
#     assert len(atoms._diffs) >= 2
    
#     # Check for correct structure in the most recent diff
#     last_diff = atoms._diffs[-1]
#     assert 'diff' in last_diff
#     assert 'timestamp' in last_diff
    
#     # Check for position changes in the diff
#     assert 'values_changed' in last_diff['diff']
    
#     # Make a more complex change with atom addition
#     atoms_copy = ATAtoms(molecule('H2O'), batch_size=100)
#     initial_diff_count = len(atoms_copy._diffs)
    
#     # Add a new atom
#     atoms_copy += Atom('H', position=(3.0, 3.0, 3.0))
    
#     # Verify diff count increased
#     assert len(atoms_copy._diffs) > initial_diff_count
    
#     # Verify diff captures added atom
#     last_diff = atoms_copy._diffs[-1]
#     assert 'iterable_item_added' in last_diff['diff'] or 'values_changed' in last_diff['diff']
    
#     # Test structure modification
#     atoms = ATAtoms(molecule('CH4'), batch_size=100)
#     initial_diff_count = len(atoms._diffs)
    
#     # Randomly displace atoms
#     atoms.rattle(stdev=0.1)
    
#     # Check if diff contains positions changes
#     assert len(atoms._diffs) > initial_diff_count
    
#     # Check consistent structure of diffs
#     for diff in atoms._diffs:
#         assert 'timestamp' in diff
#         assert ('state' in diff) != ('diff' in diff)


# def test_multiple_changes_tracking():
#     """Test tracking of multiple changes to atoms including positions, forces, etc."""
#     # Create a new system with batch_size to avoid syncing
#     atoms = ATAtoms(molecule('H2O'), batch_size=100)
    
#     # Record initial diff count
#     initial_diff_count = len(atoms._diffs)
    
#     # Perform a series of different operations
    
#     # 1. Translate the system
#     atoms.translate([1.0, 0.5, 0.2])
    
#     # 2. Set initial charges
#     atoms.set_initial_charges([0.5, -0.25, -0.25])
    
#     # 3. Apply constraint
#     from ase.constraints import FixAtoms
#     constraint = FixAtoms(indices=[0])  # Fix the first atom
#     atoms.set_constraint(constraint)
    
#     # 4. Change cell
#     atoms.set_cell([10.0, 10.0, 10.0, 90, 90, 90], scale_atoms=False)
    
#     # 5. Directly modify a position
#     atoms[2].position += [0.1, 0.1, 0.1]
    
#     # 6. Add a calculator and compute forces
#     atoms.calc = EMT()
#     forces = atoms.get_forces()
    
#     # 7. Perform a small optimization
#     opt = BFGS(atoms)
#     opt.run(steps=2)
    
#     # Now check if all these operations were tracked
#     final_diff_count = len(atoms._diffs)
#     assert final_diff_count > initial_diff_count  # Just check that we have more diffs than we started with
    
#     # Check that the diffs contain different types of changes
#     diff_types = set()
#     for diff in atoms._diffs[initial_diff_count:]:
#         if 'diff' in diff:
#             for change_type in diff['diff']:
#                 diff_types.add(change_type)
    
#     # Expect to see various types of changes
#     assert len(diff_types) >= 1  # At least one type of change


def reconstruct_atoms_from_diffs(initial_state, diffs):
    """
    Reconstruct an ATAtoms object from initial state and a series of diffs.
    
    Parameters:
    -----------
    initial_state : dict
        The initial state of the atoms object
    diffs : list
        List of diff dictionaries as stored in ATAtoms._diffs
        
    Returns:
    --------
    ATAtoms
        Reconstructed atoms object with the final state
    """
    from copy import deepcopy
    from ase import Atoms
    import numpy as np
    
    # Get the key properties from the initial state
    positions = None
    symbols = None
    pbc = None
    cell = None
    
    # Check initial state for properties
    if isinstance(initial_state, dict):
        if 'positions' in initial_state:
            positions = np.array(initial_state['positions'])
        if 'symbols' in initial_state:
            symbols = initial_state['symbols']
        if 'pbc' in initial_state:
            pbc = initial_state['pbc']
        if 'cell' in initial_state:
            cell = initial_state['cell']
    
    # Now apply all the diffs that contain full states
    for diff_entry in diffs:
        if 'state' in diff_entry and isinstance(diff_entry['state'], dict):
            state = diff_entry['state']
            if 'positions' in state:
                positions = np.array(state['positions'])
            if 'symbols' in state:
                symbols = state['symbols']
            if 'pbc' in state:
                pbc = state['pbc']
            if 'cell' in state:
                cell = state['cell']
    
    # Now apply all the value change diffs
    for diff_entry in diffs:
        if 'diff' in diff_entry and 'values_changed' in diff_entry['diff']:
            changes = diff_entry['diff']['values_changed']
            
            # Look for position changes
            for path, change in changes.items():
                if 'positions' in path and 'new_value' in change:
                    # Extract index from path like "root['positions'][0][1]"
                    parts = path.split('[')
                    if len(parts) >= 3:
                        try:
                            idx1 = int(parts[2].strip(']'))
                            if len(parts) >= 4:
                                idx2 = int(parts[3].strip(']'))
                                # Update specific position component
                                if positions is not None and idx1 < len(positions):
                                    positions[idx1][idx2] = change['new_value']
                            else:
                                # Update entire position vector
                                if positions is not None and idx1 < len(positions):
                                    positions[idx1] = change['new_value']
                        except (ValueError, IndexError):
                            continue
    
    # Create the reconstructed atoms object
    if symbols is not None and positions is not None:
        atoms = Atoms(symbols=symbols, positions=positions)
        if cell is not None:
            atoms.set_cell(cell)
        if pbc is not None:
            atoms.set_pbc(pbc)
        return ATAtoms(atoms)
    else:
        raise


# def test_reconstruct_from_diffs():
#     """Test reconstructing an atoms object from initial state and diffs."""
#     # Create a system and perform some changes
#     atoms = ATAtoms(molecule('H2O'), batch_size=100)
    
#     # Record initial state
#     initial_state = atoms._diffs[0]['state']
    
#     # Ensure we have the expected properties in the initial state
#     assert 'positions' in initial_state
#     assert 'symbols' in initial_state
    
#     # Perform some operations
#     original_positions = atoms.positions.copy()
#     atoms.translate([1.0, 0.0, 0.0])
#     assert not np.allclose(atoms.positions, original_positions)
    
#     atoms[0].position += [0.5, 0.5, 0.5]
    
#     # Get all diffs after the initial state
#     diffs = atoms._diffs[1:]
    
#     # Make sure we have some diffs
#     assert len(diffs) > 0
    
#     # Reconstruct the atoms object
#     reconstructed = reconstruct_atoms_from_diffs(initial_state, diffs)
    
#     # Verify the reconstruction
#     assert len(reconstructed) == len(atoms), f"Expected {len(atoms)} atoms, got {len(reconstructed)}"
    
#     # Check positions with a small tolerance
#     positions_match = np.allclose(reconstructed.positions, atoms.positions, atol=1e-6)
#     assert positions_match, f"\nExpected: {atoms.positions}\nGot: {reconstructed.positions}"
    
#     # Chemical formula should match exactly
#     assert reconstructed.get_chemical_formula() == atoms.get_chemical_formula(), \
#         f"Expected formula {atoms.get_chemical_formula()}, got {reconstructed.get_chemical_formula()}"


def test_save_current_state():
    """Test the save_current_state method."""
    # Create a simple atoms object
    atoms = ATAtoms(molecule('H2O'))
    
    # Make some modifications
    initial_positions = atoms.get_positions().copy()
    atoms.translate([1.0, 0.5, 0.2])
    atoms[0].position += [0.2, 0.1, 0.3]
    
    # Save the current state
    state_id = atoms.save_current_state()
    
    # Verify a state_id was returned
    assert state_id is not None
    
    # Make additional changes
    atoms.rotate(45, 'z')
    
    # Verify the positions are different from the saved state
    assert not np.allclose(atoms.get_positions(), initial_positions + [1.0, 0.5, 0.2] + [0.2, 0.1, 0.3])
    
    # Additional test could verify retrieval of the saved state 
    # if the API provides access to the saved states


def test_bfgs_optimization():
    """Test saving state during an optimization workflow."""
    # Create system and add some random displacements
    atoms = ATAtoms(bulk('Cu', cubic=True), project_id="ad7a74f8-e2b2-426c-9dc6-7471aaa19e2a")
    # hardcode for now
    atoms.rattle(stdev=0.1)  # Displace atoms randomly to create forces
    
    # Save initial positions
    initial_positions = atoms.get_positions().copy()
    
    # Setup calculator and optimizer
    atoms.calc = EMT()
    opt = BFGS(atoms)
    
    # Attach the track_changes method to the optimizer
    opt.attach(atoms.track_changes, interval=1)
    
    # Run a few steps of optimization
    opt.run(steps=2)
    
    # Save the state after initial optimization
    # state_id = atoms.save_current_state()
    
    # Check that positions changed from initial
    assert not np.allclose(atoms.get_positions(), initial_positions)
    
    # Save positions after first optimization
    intermediate_positions = atoms.get_positions().copy()
    
    # Continue optimization
    opt.run(steps=2)
    
    # Verify that positions changed after additional optimization
    assert not np.allclose(atoms.get_positions(), intermediate_positions)
    
    # Check that forces were lowered
    forces = atoms.get_forces()
    max_force = np.max(np.abs(forces))
    assert max_force < 1.0  # Using a relaxed criterion
