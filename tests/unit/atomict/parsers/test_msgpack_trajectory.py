import os
import tempfile
import pytest
import numpy as np
import time
from unittest import mock

from ase import Atoms
from ase.build import molecule
from ase.io import Trajectory as ASETrajectory
from atomict.io.trajectory import Trajectory
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def empty_trajectory(temp_dir):
    """Create an empty trajectory file path."""
    return os.path.join(temp_dir, "test_trajectory.mpk")


@pytest.fixture
def water():
    """Create a water molecule for testing."""
    return molecule('H2O')


@pytest.fixture
def water_with_calc():
    """Create a water molecule with attached calculator results."""
    water = molecule('H2O')
    
    # Create a calculator with some results
    energy = -10.0
    forces = np.array([[0.0, 0.0, 0.1], [0.0, 0.1, 0.0], [0.1, 0.0, 0.0]])
    
    calc = SinglePointCalculator(water, energy=energy, forces=forces)
    water.calc = calc
    
    return water


@pytest.fixture
def constrained_water():
    """Create a water molecule with constraints."""
    water = molecule('H2O')
    constraint = FixAtoms(indices=[0])  # Fix the oxygen atom
    water.constraints = [constraint]
    return water


@pytest.fixture
def atoms_list():
    """Create a list of atoms with same number of atoms for multi-frame testing."""
    # This ensures reshape works correctly
    h2o_1 = molecule('H2O')
    h2o_2 = molecule('H2O')
    h2o_3 = molecule('H2O')
    return [h2o_1, h2o_2, h2o_3]


# Mock for get_cell method to avoid Cell.__array__ issue if needed
def mock_get_cell(atoms_list):
    """Create a patch to replace get_cell with a version that returns a numpy array directly."""
    from ase.cell import Cell
    orig_array = Cell.__array__
    
    def patched_array(self, dtype=None, copy=True):
        return np.array(self.array, dtype=dtype, copy=copy)
        
    # Monkey patch the Cell.__array__ method
    Cell.__array__ = patched_array
    
    # Return cleanup function
    def cleanup():
        Cell.__array__ = orig_array
        
    return cleanup


@pytest.fixture(autouse=True)
def setup_cell_fix():
    """Automatically applied fixture to fix Cell.__array__ issue."""
    # Apply the cell fix if needed (depends on NumPy version)
    try:
        cleanup = mock_get_cell([molecule('H2O')])
        yield
        cleanup()
    except:
        # If patching failed, still allow tests to run with our fixed implementation
        yield


def test_trajectory_write_read_basic(water, empty_trajectory):
    """Test basic write and read operations with Trajectory."""
    # Write a single frame
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
    
    # Read it back
    with Trajectory(empty_trajectory, 'r') as traj:
        assert len(traj) == 1
        loaded_atoms = traj[0]
        
        # Check basic properties
        assert len(loaded_atoms) == len(water)
        assert loaded_atoms.get_chemical_symbols() == water.get_chemical_symbols()
        assert np.allclose(loaded_atoms.get_positions(), water.get_positions())


def test_trajectory_multiple_frames(atoms_list, empty_trajectory):
    """Test writing and reading multiple frames."""
    # Write multiple frames
    with Trajectory(empty_trajectory, 'w') as traj:
        for atoms in atoms_list:
            traj.write(atoms)
    
    # Read them back
    with Trajectory(empty_trajectory, 'r') as traj:
        assert len(traj) == len(atoms_list)
        
        # Check each frame
        for i, orig_atoms in enumerate(atoms_list):
            loaded_atoms = traj[i]
            assert len(loaded_atoms) == len(orig_atoms)
            assert loaded_atoms.get_chemical_symbols() == orig_atoms.get_chemical_symbols()
            assert np.allclose(loaded_atoms.get_positions(), orig_atoms.get_positions())


def test_trajectory_flush_interval(water, temp_dir):
    """Test that flush_interval correctly controls write frequency."""
    # First, test that the first frame is always written regardless of flush_interval
    file1 = os.path.join(temp_dir, "test_first_frame.mpk")
    traj1 = Trajectory(file1, 'w', flush_interval=100)  # Very high interval
    traj1.write(water)
    assert os.path.exists(file1), "First frame should always be written immediately"
    
    # Now test the actual flush interval behavior with a new file
    file2 = os.path.join(temp_dir, "test_interval.mpk")
    traj2 = Trajectory(file2, 'w', flush_interval=3)
    
    # First frame is always written
    traj2.write(water)
    assert os.path.exists(file2), "First frame should be written"
    
    # Delete the file to check if the next frames are written
    if os.path.exists(file2):
        os.unlink(file2)
    
    # Writing frame 2 and 3 (flush interval is 3, so we should see the file after frame 3)
    traj2.write(water)  # Frame 2
    assert not os.path.exists(file2), "Second frame shouldn't be written yet"
    
    traj2.write(water)  # Frame 3 (should trigger flush)
    assert os.path.exists(file2), "Third frame should trigger a flush"
    
    traj1.close()
    traj2.close()


def test_trajectory_append(water, atoms_list, empty_trajectory):
    """Test appending to an existing trajectory."""
    # Write initial frame
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
    
    # Append more frames
    with Trajectory(empty_trajectory, 'a') as traj:
        for atoms in atoms_list:
            traj.write(atoms)
    
    # Read and verify
    with Trajectory(empty_trajectory, 'r') as traj:
        assert len(traj) == 1 + len(atoms_list)
        
        # First frame should be water
        first = traj[0]
        assert first.get_chemical_symbols() == water.get_chemical_symbols()
        
        # Remaining frames should be from atoms_list
        for i, orig_atoms in enumerate(atoms_list):
            loaded_atoms = traj[i+1]
            assert loaded_atoms.get_chemical_symbols() == orig_atoms.get_chemical_symbols()


def test_trajectory_context_manager(water, empty_trajectory):
    """Test that context manager properly flushes data at exit."""
    # Write without explicitly closing
    with Trajectory(empty_trajectory, 'w', flush_interval=10) as traj:
        # This normally wouldn't flush with flush_interval=10
        traj.write(water)
    
    # Verify file exists and can be read
    assert os.path.exists(empty_trajectory)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        assert len(traj) == 1
        assert traj[0].get_chemical_symbols() == water.get_chemical_symbols()


def test_trajectory_calculator_data(water_with_calc, empty_trajectory):
    """Test that calculator data is preserved."""
    # Write atoms with calculator
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water_with_calc)
    
    # Read back and check calculator data
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check that it has a calculator
        assert loaded_atoms.calc is not None
        
        # Check calculator properties
        assert loaded_atoms.calc.get_potential_energy() == water_with_calc.calc.get_potential_energy()
        assert np.allclose(loaded_atoms.calc.get_forces(), water_with_calc.calc.get_forces())


def test_trajectory_constraints(constrained_water, empty_trajectory):
    """Test that constraints are preserved."""
    # Write atoms with constraints
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(constrained_water)
    
    # Read back and check constraints
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check that constraints exist
        assert len(loaded_atoms.constraints) == 1
        
        # Check that it's the right type
        assert isinstance(loaded_atoms.constraints[0], FixAtoms)
        
        # Check that the indices match
        fixed_indices = loaded_atoms.constraints[0].get_indices()
        original_indices = constrained_water.constraints[0].get_indices()
        assert np.array_equal(fixed_indices, original_indices)


def test_trajectory_slicing(atoms_list, empty_trajectory):
    """Test trajectory slicing functionality."""
    # Write multiple frames
    with Trajectory(empty_trajectory, 'w') as traj:
        for atoms in atoms_list:
            traj.write(atoms)
    
    # Read with slicing
    with Trajectory(empty_trajectory, 'r') as traj:
        # Get a slice
        sliced = traj[1:3]
        
        # Check length
        assert len(sliced) == 2
        
        # Check content
        for i, orig_idx in enumerate(range(1, 3)):
            loaded = sliced[i]
            original = atoms_list[orig_idx]
            assert loaded.get_chemical_symbols() == original.get_chemical_symbols()
            assert np.allclose(loaded.get_positions(), original.get_positions())


def test_trajectory_description(water, empty_trajectory):
    """Test setting and retrieving the description."""
    description = {"author": "Test User", "date": "2023-01-01", "system": "Water molecule"}
    
    # Write with description
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.set_description(description)
        traj.write(water)
    
    # Read and check description
    with Trajectory(empty_trajectory, 'r') as traj:
        assert traj.description is not None
        for key, value in description.items():
            assert traj.description[key] == value


def test_trajectory_manual_flush(water, empty_trajectory):
    """Test manually flushing the trajectory."""
    import os.path
    
    # Create a trajectory with high flush_interval and ensure first frame is written
    traj = Trajectory(empty_trajectory, 'w', flush_interval=1)
    
    # Write a frame and ensure it's written
    traj.write(water)
    assert os.path.exists(empty_trajectory)
    
    # Set large flush interval
    traj.flush_interval = 100
    
    # Get file size before adding more data
    size_before = os.path.getsize(empty_trajectory)
    
    # Write another frame (should not flush)
    traj.write(water)
    
    # Manually flush
    traj.flush()
    
    # File size should be larger after flush
    size_after = os.path.getsize(empty_trajectory)
    assert size_after > size_before, "File size should increase after manual flush"
    
    # Read to verify content
    with Trajectory(empty_trajectory, 'r') as read_traj:
        assert len(read_traj) == 2
    
    # Clean up
    traj.close()


def test_drop_in_replacement(water, temp_dir):
    """Test that our Trajectory is a drop-in replacement for ASE's Trajectory."""
    # Define file paths
    ase_path = os.path.join(temp_dir, "ase_traj.traj")
    our_path = os.path.join(temp_dir, "our_traj.mpk")
    
    # Write the same data with both implementations
    with ASETrajectory(ase_path, 'w') as ase_traj:
        ase_traj.write(water)
    
    with Trajectory(our_path, 'w') as our_traj:
        our_traj.write(water)
    
    # Read back with respective implementations
    with ASETrajectory(ase_path, 'r') as ase_traj:
        ase_atoms = ase_traj[0]
    
    with Trajectory(our_path, 'r') as our_traj:
        our_atoms = our_traj[0]
    
    # Compare results
    assert ase_atoms.get_chemical_symbols() == our_atoms.get_chemical_symbols()
    assert np.allclose(ase_atoms.get_positions(), our_atoms.get_positions())
    assert np.allclose(np.array(ase_atoms.get_cell()), np.array(our_atoms.get_cell()))
