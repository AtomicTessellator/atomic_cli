import pytest
from ase import Atoms


@pytest.fixture
def mock_atoms():
    """Create a simple mock Atoms object for testing"""
    return Atoms('H2O', positions=[[0, 0, 0], [0, 0, 1], [0, 1, 0]])


@pytest.fixture
def mock_atoms_list(mock_atoms):
    """Create a list of mock Atoms objects"""
    atoms1 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 1]])
    atoms2 = Atoms('CO', positions=[[0, 0, 0], [0, 0, 1.2]])
    return [atoms1, atoms2, mock_atoms]


@pytest.fixture
def cascade_atraj_file():
    """Path to the real cascade.atraj fixture file"""
    import os
    # Get the directory containing this test file
    # tests/unit/atomict/simulation/test_utils.py
    test_file = os.path.abspath(__file__)
    # Go up to tests/unit/atomict/simulation -> atomict -> unit -> tests
    tests_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(test_file))))
    fixtures_dir = os.path.join(tests_dir, 'fixtures')
    return os.path.join(fixtures_dir, 'cascade.atraj')


class TestReadAtrajWithRealData:
    """Integration tests that actually test read_atraj with real file I/O"""
    
    def test_read_atraj_from_cascade_fixture(self, cascade_atraj_file):
        """Test read_atraj with real cascade.atraj fixture file"""
        from atomict.io.formats.atraj import read_atraj
        import os
        
        # Verify the fixture file exists
        assert os.path.exists(cascade_atraj_file), f"Fixture file not found: {cascade_atraj_file}"
        
        # Read it with read_atraj
        atoms_list, metadata = read_atraj(cascade_atraj_file)
        
        # Verify we got a list back
        assert isinstance(atoms_list, list)
        assert len(atoms_list) == 1, f"Expected 1 frame, got {len(atoms_list)}"
        
        # Verify the structure
        atoms = atoms_list[0]
        assert len(atoms) == 22, f"Expected 22 atoms, got {len(atoms)}"
        
        # Check chemical composition (16 C + 6 Ne = 22 atoms based on the data)
        formula = atoms.get_chemical_formula()
        assert 'C' in formula, f"Expected Carbon in formula, got {formula}"
        assert 'Ne' in formula, f"Expected Neon in formula, got {formula}"
        
        # Verify metadata is a dict
        assert isinstance(metadata, dict)
        
        # Verify positions array exists and has correct shape
        assert atoms.positions is not None
        assert atoms.positions.shape == (22, 3), f"Expected (22, 3) positions, got {atoms.positions.shape}"
        
        # Verify cell is present
        assert atoms.cell is not None
    
    def test_read_atraj_single_atoms(self, mock_atoms, tmp_path):
        """Test read_atraj with a real .atraj file containing a single Atoms object"""
        from atomict.io.formats.atraj import write_atraj, read_atraj
        
        # Create a test file
        atraj_file = tmp_path / "test_single.atraj"
        metadata = {'test_key': 'test_value', 'energy': -10.5}
        
        # Write the atoms object
        write_atraj(mock_atoms, str(atraj_file), metadata=metadata)
        
        # Read it back
        atoms_list, read_metadata = read_atraj(str(atraj_file))
        
        # Verify we got a list back
        assert isinstance(atoms_list, list)
        assert len(atoms_list) == 1
        
        # Verify the atoms object matches
        result_atoms = atoms_list[0]
        assert len(result_atoms) == len(mock_atoms)
        assert result_atoms.get_chemical_formula() == mock_atoms.get_chemical_formula()
        assert (result_atoms.positions == mock_atoms.positions).all()
        
        # Verify metadata
        assert read_metadata == metadata
    
    def test_read_atraj_trajectory(self, mock_atoms_list, tmp_path):
        """Test read_atraj with a real .atraj file containing multiple frames"""
        from atomict.io.formats.atraj import write_atraj, read_atraj
        
        # Create a test file with trajectory
        atraj_file = tmp_path / "test_trajectory.atraj"
        metadata = {'frames': len(mock_atoms_list), 'source': 'test'}
        
        # Write the trajectory
        write_atraj(mock_atoms_list, str(atraj_file), metadata=metadata)
        
        # Read it back
        atoms_list, read_metadata = read_atraj(str(atraj_file))
        
        # Verify we got all frames
        assert isinstance(atoms_list, list)
        assert len(atoms_list) == len(mock_atoms_list)
        
        # Verify each frame matches
        for i, (result, expected) in enumerate(zip(atoms_list, mock_atoms_list)):
            assert len(result) == len(expected), f"Frame {i} has wrong number of atoms"
            assert result.get_chemical_formula() == expected.get_chemical_formula()
        
        # Verify metadata
        assert read_metadata == metadata
    
    def test_read_atraj_with_calculated_properties(self, tmp_path):
        """Test read_atraj preserves calculated properties like energy, forces"""
        from atomict.io.formats.atraj import write_atraj, read_atraj
        from ase import Atoms
        from ase.calculators.singlepoint import SinglePointCalculator
        import numpy as np
        
        # Create atoms with calculator results
        atoms = Atoms('CO', positions=[[0, 0, 0], [0, 0, 1.2]])
        energy = -15.3
        forces = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]])
        stress = np.array([1.0, 2.0, 3.0, 0.0, 0.0, 0.0])
        
        calc = SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress)
        atoms.calc = calc
        
        # Write to file
        atraj_file = tmp_path / "test_with_calc.atraj"
        write_atraj(atoms, str(atraj_file))
        
        # Read back
        atoms_list, _ = read_atraj(str(atraj_file))
        result = atoms_list[0]
        
        # Verify calculated properties are preserved
        assert result.calc is not None
        assert abs(result.get_potential_energy() - energy) < 1e-10
        assert np.allclose(result.get_forces(), forces)
        assert np.allclose(result.get_stress(), stress)
