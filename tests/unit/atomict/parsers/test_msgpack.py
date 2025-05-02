import os
import tempfile
import pytest
import numpy as np
from unittest import mock

from ase import Atoms
from ase.build import molecule
from atomict.io.msgpack import load_msgpack, save_msgpack


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
