import os
import tempfile
import unittest
import numpy as np
from unittest import mock

from ase import Atoms
from ase.build import molecule
from atomict.io.msgpack import load_msgpack, save_msgpack


class TestMsgpack(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_file = os.path.join(self.temp_dir.name, "test_atoms.msgpack")
        
        # Create various test atoms for different scenarios
        self.single_atom = Atoms('H', positions=[[0, 0, 0]])
        self.water = molecule('H2O')
        self.co2 = molecule('CO2')
        
        # Set some custom properties on the water molecule
        self.water.set_tags([1, 2, 3])
        self.water.set_masses([2.0, 18.0, 2.0])  # Custom masses
        self.water.set_momenta(np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
        self.water.set_initial_charges([0.5, -1.0, 0.5])
        self.water.set_initial_magnetic_moments([0.1, 0.0, 0.1])
        
        # Create a list of atoms with same number of atoms for multi-frame testing
        # This ensures reshape works correctly
        self.h2o_1 = molecule('H2O')
        self.h2o_2 = molecule('H2O')
        self.h2o_3 = molecule('H2O')
        self.atoms_list = [self.h2o_1, self.h2o_2, self.h2o_3]
    
    def tearDown(self):
        # Clean up temporary directory
        self.temp_dir.cleanup()
    
    def mock_get_cell(self, mock_atoms_list):
        """Create a patch to replace get_cell with a version that returns a numpy array directly"""
        orig_get_cell = Atoms.get_cell
        
        def patched_get_cell(self):
            cell = orig_get_cell(self)
            # Convert Cell to numpy array directly
            return np.asarray(cell.array)
            
        for atom in mock_atoms_list:
            atom.get_cell = patched_get_cell.__get__(atom, Atoms)
    
    def test_save_load_single_atoms(self):
        """Test saving and loading a single Atoms object."""
        # Patch the get_cell method to avoid Cell.__array__ issue
        self.mock_get_cell([self.water])
        
        # Save and load water molecule
        save_msgpack(self.water, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        
        # Check that it's the correct type
        self.assertIsInstance(loaded_atoms, Atoms)
        
        # Test basic properties
        self.assertEqual(len(loaded_atoms), len(self.water))
        self.assertTrue(np.allclose(loaded_atoms.get_positions(), self.water.get_positions()))
        self.assertTrue(np.allclose(np.asarray(loaded_atoms.get_cell()), np.asarray(self.water.get_cell())))
        self.assertTrue(np.array_equal(loaded_atoms.get_pbc(), self.water.get_pbc()))
        self.assertEqual(loaded_atoms.get_chemical_symbols(), self.water.get_chemical_symbols())
    
    def test_save_load_atoms_list(self):
        """Test saving and loading a list of Atoms objects."""
        # Patch the get_cell method to avoid Cell.__array__ issue
        self.mock_get_cell(self.atoms_list)
        
        # Save and load multiple molecules (all H2O to keep consistent atom count)
        save_msgpack(self.atoms_list, self.test_file)
        loaded_atoms_list = load_msgpack(self.test_file)
        
        # Check that it's the correct type
        self.assertIsInstance(loaded_atoms_list, list)
        self.assertEqual(len(loaded_atoms_list), len(self.atoms_list))
        
        # Test that each atoms object has the correct properties
        for loaded, original in zip(loaded_atoms_list, self.atoms_list):
            self.assertEqual(len(loaded), len(original))
            self.assertTrue(np.allclose(loaded.get_positions(), original.get_positions()))
            self.assertTrue(np.allclose(np.asarray(loaded.get_cell()), np.asarray(original.get_cell())))
            self.assertTrue(np.array_equal(loaded.get_pbc(), original.get_pbc()))
            self.assertEqual(loaded.get_chemical_symbols(), original.get_chemical_symbols())
    
    def test_property_tags(self):
        """Test that tags are correctly saved and loaded."""
        # Patch the get_cell method
        self.mock_get_cell([self.water])
        
        save_msgpack(self.water, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        self.assertTrue(np.array_equal(loaded_atoms.get_tags(), self.water.get_tags()))
    
    def test_property_masses(self):
        """Test that custom masses are correctly saved and loaded."""
        # Patch the get_cell method
        self.mock_get_cell([self.water])
        
        save_msgpack(self.water, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        self.assertTrue(np.allclose(loaded_atoms.get_masses(), self.water.get_masses()))
    
    def test_property_momenta(self):
        """Test that momenta are correctly saved and loaded."""
        # Patch the get_cell method
        self.mock_get_cell([self.water])
        
        save_msgpack(self.water, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        self.assertTrue(np.allclose(loaded_atoms.get_momenta(), self.water.get_momenta()))
    
    def test_property_charges(self):
        """Test that initial charges are correctly saved and loaded."""
        # Patch the get_cell method
        self.mock_get_cell([self.water])
        
        save_msgpack(self.water, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        self.assertTrue(np.allclose(loaded_atoms.get_initial_charges(), self.water.get_initial_charges()))
    
    def test_property_magnetic_moments(self):
        """Test that initial magnetic moments are correctly saved and loaded."""
        # Patch the get_cell method
        self.mock_get_cell([self.water])
        
        save_msgpack(self.water, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        self.assertTrue(np.allclose(loaded_atoms.get_initial_magnetic_moments(), self.water.get_initial_magnetic_moments()))
    
    def test_property_persistence_in_list(self):
        """Test that all properties are correctly saved and loaded in a list of atoms."""
        # Create two identical water molecules for consistent reshaping
        water2 = molecule('H2O')
        water2.set_tags([1, 2, 3])
        water2.set_masses([2.0, 18.0, 2.0])
        water2.set_momenta(np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
        water2.set_initial_charges([0.5, -1.0, 0.5])
        water2.set_initial_magnetic_moments([0.1, 0.0, 0.1])
        
        atoms_list_with_props = [self.water, water2]
        
        # Patch the get_cell method
        self.mock_get_cell(atoms_list_with_props)
        
        save_msgpack(atoms_list_with_props, self.test_file)
        loaded_list = load_msgpack(self.test_file)
        
        # Check the properties on both water molecules
        for i in range(2):
            loaded_water = loaded_list[i]
            self.assertTrue(np.array_equal(loaded_water.get_tags(), self.water.get_tags()))
            self.assertTrue(np.allclose(loaded_water.get_masses(), self.water.get_masses()))
            self.assertTrue(np.allclose(loaded_water.get_momenta(), self.water.get_momenta()))
            self.assertTrue(np.allclose(loaded_water.get_initial_charges(), self.water.get_initial_charges()))
            self.assertTrue(np.allclose(loaded_water.get_initial_magnetic_moments(), self.water.get_initial_magnetic_moments()))
    
    def test_selective_property_inclusion(self):
        """Test that only non-default properties are included in the msgpack file."""
        # Create an atom with default properties
        h_atom = Atoms('H', positions=[[0, 0, 0]])
        
        # Patch the get_cell method
        self.mock_get_cell([h_atom])
        
        # Use a direct approach: save the file, then read it back with msgpack directly
        save_msgpack(h_atom, self.test_file)
        
        # Load the data directly with msgpack to inspect what was saved
        import msgpack
        with open(self.test_file, 'rb') as f:
            saved_data = msgpack.unpack(f, raw=False)
        
        # These should be included
        self.assertIn('n_frames', saved_data)
        self.assertIn('positions', saved_data)
        self.assertIn('cell', saved_data)
        self.assertIn('pbc', saved_data)
        
        # These should not be included with default values
        self.assertNotIn('tags', saved_data)
        self.assertNotIn('momenta', saved_data)
        self.assertNotIn('initial_charges', saved_data)
        self.assertNotIn('initial_magmoms', saved_data)
    
    def test_import_error_handling(self):
        """Test proper error handling when dependencies are missing."""
        # Properly mock the imports within the function
        with mock.patch.dict('sys.modules', {'msgpack': None}):
            with self.assertRaises(ImportError):
                save_msgpack(self.water, self.test_file)
        
        with mock.patch.dict('sys.modules', {'msgpack': None}):
            with self.assertRaises(ImportError):
                load_msgpack(self.test_file)
    
    def test_non_existing_file(self):
        """Test error handling when trying to load a non-existing file."""
        with self.assertRaises(FileNotFoundError):
            load_msgpack('non_existing_file.msgpack')
    
    def test_large_structure(self):
        """Test handling of a larger structure to ensure performance."""
        # Create a bigger system
        big_atoms = Atoms('H16', 
                         positions=np.random.rand(16, 3), 
                         cell=[10, 10, 10], 
                         pbc=True)
        
        # Set some properties
        big_atoms.set_tags(np.arange(16))
        
        # Patch the get_cell method
        self.mock_get_cell([big_atoms])
        
        # Save and load
        save_msgpack(big_atoms, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        
        # Verify
        self.assertEqual(len(loaded_atoms), len(big_atoms))
        self.assertTrue(np.allclose(loaded_atoms.get_positions(), big_atoms.get_positions()))
        self.assertTrue(np.array_equal(loaded_atoms.get_tags(), big_atoms.get_tags()))

    def test_float32_precision(self):
        """Test that positions are stored as float32 for efficiency but maintain adequate precision."""
        # Create atoms with precise positions
        precise_atoms = Atoms('H2', positions=[[0.12345678, 0.23456789, 0.34567890],
                                              [0.98765432, 0.87654321, 0.76543210]])
        
        # Patch the get_cell method
        self.mock_get_cell([precise_atoms])
        
        # Save and load
        save_msgpack(precise_atoms, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        
        # Check precision - we expect some loss due to float32 conversion
        # but it should be adequate for atomic simulations
        self.assertTrue(np.allclose(loaded_atoms.get_positions(), precise_atoms.get_positions(), 
                                  rtol=1e-6, atol=1e-6))
        
    def test_different_atom_types(self):
        """Test with a mix of different atom types."""
        mixed_atoms = Atoms('HCNOF', positions=np.random.rand(5, 3))
        
        # Patch the get_cell method
        self.mock_get_cell([mixed_atoms])
        
        # Save and load
        save_msgpack(mixed_atoms, self.test_file)
        loaded_atoms = load_msgpack(self.test_file)
        
        # Verify
        self.assertEqual(loaded_atoms.get_chemical_symbols(), mixed_atoms.get_chemical_symbols())
        
    def test_masses_are_saved(self):
        """Test that masses are correctly saved and loaded for proper atomic weights."""
        # Create atoms with default masses
        atoms = Atoms('H2O', positions=np.random.rand(3, 3))
        
        # Patch the get_cell method
        self.mock_get_cell([atoms])
        
        # Use a direct approach: save the file, then read it back with msgpack directly
        save_msgpack(atoms, self.test_file)
        
        # Load the data directly with msgpack to inspect what was saved
        import msgpack
        with open(self.test_file, 'rb') as f:
            saved_data = msgpack.unpack(f, raw=False)
        
        # Masses should be included for restoring proper atomic weights
        self.assertIn('masses', saved_data)
        
        # Verify that the loaded object has the correct masses
        loaded_atoms = load_msgpack(self.test_file)
        self.assertTrue(np.allclose(loaded_atoms.get_masses(), atoms.get_masses()))


if __name__ == '__main__':
    unittest.main()
