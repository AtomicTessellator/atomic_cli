import os
import tempfile
import pytest
import numpy as np
from unittest import mock

from ase import Atoms
from ase.build import molecule
from atomict.io.msgpack import load_msgpack, save_msgpack
from ase.io import Trajectory
from ase.constraints import FixAtoms

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def test_file(temp_dir):
    """Create a test file path."""
    return os.path.join(temp_dir, "test_atoms.msgpack")


@pytest.fixture
def single_atom():
    """Create a single atom for testing."""
    return Atoms('H', positions=[[0, 0, 0]])


@pytest.fixture
def water():
    """Create a water molecule with custom properties."""
    water = molecule('H2O')
    water.set_tags([1, 2, 3])
    water.set_masses([2.0, 18.0, 2.0])  # Custom masses
    water.set_momenta(np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    water.set_initial_charges([0.5, -1.0, 0.5])
    water.set_initial_magnetic_moments([0.1, 0.0, 0.1])
    return water


@pytest.fixture
def co2():
    """Create a CO2 molecule for testing."""
    return molecule('CO2')


@pytest.fixture
def atoms_list():
    """Create a list of atoms with same number of atoms for multi-frame testing."""
    # This ensures reshape works correctly
    h2o_1 = molecule('H2O')
    h2o_2 = molecule('H2O')
    h2o_3 = molecule('H2O')
    return [h2o_1, h2o_2, h2o_3]


@pytest.fixture
def atoms_with_all_props():
    """Create an Atoms object with all the properties we want to test."""
    atoms = molecule('H2O')
    
    # Set standard properties
    atoms.set_tags([1, 2, 3])
    atoms.set_masses([2.0, 18.0, 2.0])
    atoms.set_momenta(np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    atoms.set_initial_charges([0.5, -1.0, 0.5])
    atoms.set_initial_magnetic_moments([0.1, 0.0, 0.1])
    
    # Set custom attributes
    atoms.ase_objtype = 'atoms'
    atoms.top_mask = np.array([False, True, False], dtype=bool)
    
    # Don't set number_of_lattice_vectors directly since it's a property without a setter
    # Instead, adjust the cell which affects this property
    atoms.set_cell([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]])
    atoms.set_pbc(True)
    
    # Add constraint
    constraint = FixAtoms(indices=[0])
    atoms.constraints = [constraint]
    
    # Set forces and stress
    forces = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    atoms.arrays['forces'] = forces
    atoms.stress = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    
    return atoms


def mock_get_cell(atoms_list):
    """Create a patch to replace get_cell with a version that returns a numpy array directly."""
    orig_get_cell = Atoms.get_cell
    
    def patched_get_cell(self):
        cell = orig_get_cell(self)
        # Convert Cell to numpy array directly
        return np.asarray(cell.array)
        
    for atom in atoms_list:
        atom.get_cell = patched_get_cell.__get__(atom, Atoms)


def test_save_load_single_atoms(water, test_file):
    """Test saving and loading a single Atoms object."""
    # Patch the get_cell method to avoid Cell.__array__ issue
    mock_get_cell([water])
    
    # Save and load water molecule
    save_msgpack(water, test_file)
    loaded_atoms = load_msgpack(test_file)
    
    # Check that it's the correct type
    assert isinstance(loaded_atoms, Atoms)
    
    # Test basic properties
    assert len(loaded_atoms) == len(water)
    assert np.allclose(loaded_atoms.get_positions(), water.get_positions())
    assert np.allclose(np.asarray(loaded_atoms.get_cell()), np.asarray(water.get_cell()))
    assert np.array_equal(loaded_atoms.get_pbc(), water.get_pbc())
    assert loaded_atoms.get_chemical_symbols() == water.get_chemical_symbols()


def test_save_load_atoms_list(atoms_list, test_file):
    """Test saving and loading a list of Atoms objects."""
    # Patch the get_cell method to avoid Cell.__array__ issue
    mock_get_cell(atoms_list)
    
    # Save and load multiple molecules (all H2O to keep consistent atom count)
    save_msgpack(atoms_list, test_file)
    loaded_atoms_list = load_msgpack(test_file)
    
    # Check that it's the correct type
    assert isinstance(loaded_atoms_list, list)
    assert len(loaded_atoms_list) == len(atoms_list)
    
    # Test that each atoms object has the correct properties
    for loaded, original in zip(loaded_atoms_list, atoms_list):
        assert len(loaded) == len(original)
        assert np.allclose(loaded.get_positions(), original.get_positions())
        assert np.allclose(np.asarray(loaded.get_cell()), np.asarray(original.get_cell()))
        assert np.array_equal(loaded.get_pbc(), original.get_pbc())
        assert loaded.get_chemical_symbols() == original.get_chemical_symbols()


def test_property_tags(water, test_file):
    """Test that tags are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([water])
    
    save_msgpack(water, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert np.array_equal(loaded_atoms.get_tags(), water.get_tags())


def test_property_masses(water, test_file):
    """Test that custom masses are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([water])
    
    save_msgpack(water, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert np.allclose(loaded_atoms.get_masses(), water.get_masses())


def test_property_momenta(water, test_file):
    """Test that momenta are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([water])
    
    save_msgpack(water, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert np.allclose(loaded_atoms.get_momenta(), water.get_momenta())


def test_property_charges(water, test_file):
    """Test that initial charges are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([water])
    
    save_msgpack(water, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert np.allclose(loaded_atoms.get_initial_charges(), water.get_initial_charges())


def test_property_magnetic_moments(water, test_file):
    """Test that initial magnetic moments are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([water])
    
    save_msgpack(water, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert np.allclose(loaded_atoms.get_initial_magnetic_moments(), water.get_initial_magnetic_moments())


def test_property_persistence_in_list(water, test_file):
    """Test that all properties are correctly saved and loaded in a list of atoms."""
    # Create two identical water molecules for consistent reshaping
    water2 = molecule('H2O')
    water2.set_tags([1, 2, 3])
    water2.set_masses([2.0, 18.0, 2.0])
    water2.set_momenta(np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    water2.set_initial_charges([0.5, -1.0, 0.5])
    water2.set_initial_magnetic_moments([0.1, 0.0, 0.1])
    
    atoms_list_with_props = [water, water2]
    
    # Patch the get_cell method
    mock_get_cell(atoms_list_with_props)
    
    save_msgpack(atoms_list_with_props, test_file)
    loaded_list = load_msgpack(test_file)
    
    # Check the properties on both water molecules
    for i in range(2):
        loaded_water = loaded_list[i]
        assert np.array_equal(loaded_water.get_tags(), water.get_tags())
        assert np.allclose(loaded_water.get_masses(), water.get_masses())
        assert np.allclose(loaded_water.get_momenta(), water.get_momenta())
        assert np.allclose(loaded_water.get_initial_charges(), water.get_initial_charges())
        assert np.allclose(loaded_water.get_initial_magnetic_moments(), water.get_initial_magnetic_moments())


def test_selective_property_inclusion(test_file):
    """Test that only non-default properties are included in the msgpack file."""
    # Create an atom with default properties
    h_atom = Atoms('H', positions=[[0, 0, 0]])
    
    # Patch the get_cell method
    mock_get_cell([h_atom])
    
    # Use a direct approach: save the file, then read it back with msgpack directly
    save_msgpack(h_atom, test_file)
    
    # Load the data directly with msgpack to inspect what was saved
    import msgpack
    with open(test_file, 'rb') as f:
        saved_data = msgpack.unpack(f, raw=False)
    
    # These should be included
    assert 'n_frames' in saved_data
    assert 'positions' in saved_data
    assert 'cell' in saved_data
    assert 'pbc' in saved_data
    
    # These should not be included with default values
    assert 'tags' not in saved_data
    assert 'momenta' not in saved_data
    assert 'initial_charges' not in saved_data
    assert 'initial_magmoms' not in saved_data


def test_import_error_handling():
    """Test proper error handling when dependencies are missing."""
    # Properly mock the imports within the function
    with mock.patch.dict('sys.modules', {'msgpack': None}):
        with pytest.raises(ImportError):
            save_msgpack(molecule('H2O'), 'dummy_file.msgpack')
    
    with mock.patch.dict('sys.modules', {'msgpack': None}):
        with pytest.raises(ImportError):
            load_msgpack('dummy_file.msgpack')


def test_non_existing_file():
    """Test error handling when trying to load a non-existing file."""
    with pytest.raises(FileNotFoundError):
        load_msgpack('non_existing_file.msgpack')


def test_large_structure(test_file):
    """Test handling of a larger structure to ensure performance."""
    # Create a bigger system
    big_atoms = Atoms('H16', 
                     positions=np.random.rand(16, 3), 
                     cell=[10, 10, 10], 
                     pbc=True)
    
    # Set some properties
    big_atoms.set_tags(np.arange(16))
    
    # Patch the get_cell method
    mock_get_cell([big_atoms])
    
    # Save and load
    save_msgpack(big_atoms, test_file)
    loaded_atoms = load_msgpack(test_file)
    
    # Verify
    assert len(loaded_atoms) == len(big_atoms)
    assert np.allclose(loaded_atoms.get_positions(), big_atoms.get_positions())
    assert np.array_equal(loaded_atoms.get_tags(), big_atoms.get_tags())


def test_float32_precision(test_file):
    """Test that positions are stored as float32 for efficiency but maintain adequate precision."""
    # Create atoms with precise positions
    precise_atoms = Atoms('H2', positions=[[0.12345678, 0.23456789, 0.34567890],
                                          [0.98765432, 0.87654321, 0.76543210]])
    
    # Patch the get_cell method
    mock_get_cell([precise_atoms])
    
    # Save and load
    save_msgpack(precise_atoms, test_file)
    loaded_atoms = load_msgpack(test_file)
    
    # Check precision - we expect some loss due to float32 conversion
    # but it should be adequate for atomic simulations
    assert np.allclose(loaded_atoms.get_positions(), precise_atoms.get_positions(), 
                      rtol=1e-6, atol=1e-6)
    
    
def test_different_atom_types(test_file):
    """Test with a mix of different atom types."""
    mixed_atoms = Atoms('HCNOF', positions=np.random.rand(5, 3))
    
    # Patch the get_cell method
    mock_get_cell([mixed_atoms])
    
    # Save and load
    save_msgpack(mixed_atoms, test_file)
    loaded_atoms = load_msgpack(test_file)
    
    # Verify
    assert loaded_atoms.get_chemical_symbols() == mixed_atoms.get_chemical_symbols()
    
    
def test_masses_are_saved(test_file):
    """Test that masses are correctly saved and loaded for proper atomic weights."""
    # Create atoms with default masses
    atoms = Atoms('H2O', positions=np.random.rand(3, 3))
    
    # Patch the get_cell method
    mock_get_cell([atoms])
    
    # Use a direct approach: save the file, then read it back with msgpack directly
    save_msgpack(atoms, test_file)
    
    # Load the data directly with msgpack to inspect what was saved
    import msgpack
    with open(test_file, 'rb') as f:
        saved_data = msgpack.unpack(f, raw=False)
    
    # Masses should be included for restoring proper atomic weights
    assert 'masses' in saved_data
    
    # Verify that the loaded object has the correct masses
    loaded_atoms = load_msgpack(test_file)
    assert np.allclose(loaded_atoms.get_masses(), atoms.get_masses())


def test_property_ase_objtype(atoms_with_all_props, test_file):
    """Test that ase_objtype is correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert hasattr(loaded_atoms, 'ase_objtype')
    assert loaded_atoms.ase_objtype == atoms_with_all_props.ase_objtype


def test_property_top_mask(atoms_with_all_props, test_file):
    """Test that top_mask is correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert hasattr(loaded_atoms, 'top_mask')
    assert np.array_equal(loaded_atoms.top_mask, atoms_with_all_props.top_mask)


def test_property_number_of_lattice_vectors(atoms_with_all_props, test_file):
    """Test that number_of_lattice_vectors is correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    
    # Check cell rank instead of deprecated number_of_lattice_vectors
    assert loaded_atoms.cell.rank == atoms_with_all_props.cell.rank


def test_property_constraints(atoms_with_all_props, test_file):
    """Test that constraints are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    
    assert len(loaded_atoms.constraints) == len(atoms_with_all_props.constraints)
    assert isinstance(loaded_atoms.constraints[0], type(atoms_with_all_props.constraints[0]))
    
    # Check constraint details
    original_constraint = atoms_with_all_props.constraints[0]
    loaded_constraint = loaded_atoms.constraints[0]
    
    if isinstance(original_constraint, FixAtoms):
        assert np.array_equal(original_constraint.index, loaded_constraint.index)


def test_property_numbers(atoms_with_all_props, test_file):
    """Test that atomic numbers are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert np.array_equal(loaded_atoms.get_atomic_numbers(), atoms_with_all_props.get_atomic_numbers())


def test_property_forces(atoms_with_all_props, test_file):
    """Test that forces are correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert 'forces' in loaded_atoms.arrays
    assert np.allclose(loaded_atoms.arrays['forces'], atoms_with_all_props.arrays['forces'])


def test_property_stress(atoms_with_all_props, test_file):
    """Test that stress is correctly saved and loaded."""
    # Patch the get_cell method
    mock_get_cell([atoms_with_all_props])
    
    save_msgpack(atoms_with_all_props, test_file)
    loaded_atoms = load_msgpack(test_file)
    assert hasattr(loaded_atoms, 'stress')
    assert np.allclose(loaded_atoms.stress, atoms_with_all_props.stress)


def test_all_properties_list(atoms_with_all_props, test_file):
    """Test that all properties are correctly saved and loaded in a list context."""
    # Create a second copy of the atoms with all properties
    atoms2 = atoms_with_all_props.copy()
    # Ensure custom attributes are copied correctly
    atoms2.ase_objtype = atoms_with_all_props.ase_objtype
    atoms2.top_mask = atoms_with_all_props.top_mask.copy()  # Explicitly copy the top_mask
    # Ensure stress is copied correctly
    if hasattr(atoms_with_all_props, 'stress'):
        atoms2.stress = atoms_with_all_props.stress.copy()
    
    atoms_list = [atoms_with_all_props, atoms2]
    
    # Patch the get_cell method
    mock_get_cell(atoms_list)
    
    save_msgpack(atoms_list, test_file)
    loaded_list = load_msgpack(test_file)
    
    for i in range(2):
        loaded_atoms = loaded_list[i]
        original_atoms = atoms_list[i]
        
        # Check all standard properties
        assert np.array_equal(loaded_atoms.get_tags(), original_atoms.get_tags())
        assert np.allclose(loaded_atoms.get_masses(), original_atoms.get_masses())
        assert np.allclose(loaded_atoms.get_momenta(), original_atoms.get_momenta())
        assert np.allclose(loaded_atoms.get_initial_charges(), original_atoms.get_initial_charges())
        assert np.allclose(loaded_atoms.get_initial_magnetic_moments(), original_atoms.get_initial_magnetic_moments())
        
        # Check new properties
        assert hasattr(loaded_atoms, 'ase_objtype')
        assert loaded_atoms.ase_objtype == original_atoms.ase_objtype
        
        assert hasattr(loaded_atoms, 'top_mask')
        assert np.array_equal(loaded_atoms.top_mask, original_atoms.top_mask)
        
        # Check number_of_lattice_vectors indirectly
        assert loaded_atoms.cell.rank == original_atoms.cell.rank
        
        assert len(loaded_atoms.constraints) == len(original_atoms.constraints)
        assert np.array_equal(loaded_atoms.get_atomic_numbers(), original_atoms.get_atomic_numbers())
        
        # Check forces and stress
        assert 'forces' in loaded_atoms.arrays
        assert np.allclose(loaded_atoms.arrays['forces'], original_atoms.arrays['forces'])
        
        # Ensure stress is handled correctly
        assert hasattr(loaded_atoms, 'stress')
        assert hasattr(original_atoms, 'stress')
        assert np.allclose(loaded_atoms.stress, original_atoms.stress)
