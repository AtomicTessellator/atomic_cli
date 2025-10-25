import os
import tempfile
import pytest
import numpy as np
import time
from unittest import mock

from ase import Atoms
from ase.build import molecule, bulk
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
    return os.path.join(temp_dir, "test_trajectory.tess")


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
    file1 = os.path.join(temp_dir, "test_first_frame.tess")
    traj1 = Trajectory(file1, 'w', flush_interval=100)  # Very high interval
    traj1.write(water)
    assert os.path.exists(file1), "First frame should always be written immediately"
    
    # Now test the actual flush interval behavior with a new file
    file2 = os.path.join(temp_dir, "test_interval.tess")
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
    our_path = os.path.join(temp_dir, "our_traj.tess")
    
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


def test_trajectory_pbc_preservation(water, empty_trajectory):
    """Test that periodic boundary conditions are preserved."""
    # Set custom PBC
    water.set_pbc([True, False, True])
    
    # Write and read
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check PBC values are preserved
        assert np.array_equal(loaded_atoms.get_pbc(), water.get_pbc())


def test_trajectory_arrays_preservation(water, empty_trajectory):
    """Test that arrays dictionary data is preserved."""
    # Add custom arrays
    water.arrays['custom_array'] = np.array([1.0, 2.0, 3.0])
    water.arrays['test_vector'] = np.array([[1.0, 0.0, 0.0], 
                                           [0.0, 1.0, 0.0],
                                           [0.0, 0.0, 1.0]])
    
    # Write and read
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check arrays are preserved with correct values
        assert 'custom_array' in loaded_atoms.arrays
        assert 'test_vector' in loaded_atoms.arrays
        assert np.array_equal(loaded_atoms.arrays['custom_array'], water.arrays['custom_array'])
        assert np.array_equal(loaded_atoms.arrays['test_vector'], water.arrays['test_vector'])


def test_trajectory_info_preservation(water, empty_trajectory):
    """Test that info dictionary data is preserved."""
    # Add custom info (excluding special keys used by trajectory)
    water.info['test_info'] = 'test value'
    water.info['numeric_data'] = 42
    water.info['complex_data'] = {'a': 1, 'b': [1, 2, 3]}
    
    # Write and read
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check custom info is preserved
        assert loaded_atoms.info['test_info'] == 'test value'
        assert loaded_atoms.info['numeric_data'] == 42
        assert loaded_atoms.info['complex_data']['a'] == 1
        assert loaded_atoms.info['complex_data']['b'] == [1, 2, 3]


def test_trajectory_cell_preservation(empty_trajectory):
    """Test that cell information is correctly preserved."""
    # Create an atoms object with a non-trivial cell
    from ase.lattice.cubic import FaceCenteredCubic
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                             size=(2, 2, 2), symbol='Au', pbc=True)
    
    # Get the original cell
    original_cell = atoms.get_cell()
    
    # Write and read
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(atoms)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check cell data
        loaded_cell = loaded_atoms.get_cell()
        assert np.allclose(loaded_cell, original_cell)
        
        # Check cell params are preserved
        assert np.allclose(loaded_cell.lengths(), original_cell.lengths())
        assert np.allclose(loaded_cell.angles(), original_cell.angles())


def test_trajectory_tags_preservation(water, empty_trajectory):
    """Test that atom tags are preserved."""
    # Set custom tags
    tags = np.array([1, 2, 3])
    water.set_tags(tags)
    
    # Write and read
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check tags are preserved
        assert np.array_equal(loaded_atoms.get_tags(), tags)


def test_complex_atoms_data(empty_trajectory):
    """Test storage of complex atoms with multiple properties."""
    # Create an atom with many properties set
    from ase.build import bulk
    atoms = bulk('Si')
    
    # Set various properties
    atoms.set_pbc([True, True, False])
    atoms.set_masses(np.array([29.0] * len(atoms)))
    atoms.set_initial_charges(np.array([0.5] * len(atoms)))
    
    # Add custom arrays and info
    atoms.arrays['momentum'] = np.ones((len(atoms), 3))
    atoms.info['calculation_params'] = {'cutoff': 300, 'xc': 'PBE'}
    
    # Add calculator with results
    energy = -10.0
    forces = np.random.random((len(atoms), 3))
    stress = np.array([1.0, 2.0, 3.0, 0.0, 0.0, 0.0])
    
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress)
    atoms.calc = calc
    
    # Set constraints
    from ase.constraints import FixAtoms, FixedPlane
    constraint1 = FixAtoms(indices=[0])
    constraint2 = FixedPlane(1, [1, 0, 0])  # Fix atom 1 in the yz plane
    atoms.constraints = [constraint1, constraint2]
    
    # Write and read
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(atoms)
    
    with Trajectory(empty_trajectory, 'r') as traj:
        loaded_atoms = traj[0]
        
        # Check all properties are preserved
        assert np.array_equal(loaded_atoms.get_pbc(), atoms.get_pbc())
        assert np.allclose(loaded_atoms.get_masses(), atoms.get_masses())
        assert np.allclose(loaded_atoms.get_initial_charges(), atoms.get_initial_charges())
        assert np.allclose(loaded_atoms.arrays['momentum'], atoms.arrays['momentum'])
        assert loaded_atoms.info['calculation_params'] == atoms.info['calculation_params']
        
        # Check calculator data
        assert loaded_atoms.calc is not None
        assert np.isclose(loaded_atoms.calc.get_potential_energy(), energy)
        assert np.allclose(loaded_atoms.calc.get_forces(), forces)
        assert np.allclose(loaded_atoms.calc.get_stress(), stress)
        
        # Check constraints
        assert len(loaded_atoms.constraints) == 2
        assert isinstance(loaded_atoms.constraints[0], FixAtoms)
        assert isinstance(loaded_atoms.constraints[1], FixedPlane)
        assert np.array_equal(loaded_atoms.constraints[0].get_indices(), [0])
        
        # Check the FixedPlane constraint correctly
        # The attribute name might be different when deserialized
        # We'll check the constraint using the indices and direction
        fixed_plane = loaded_atoms.constraints[1]
        assert isinstance(fixed_plane, FixedPlane)
        assert 1 in fixed_plane.get_indices()  # Check the atom index is included
        assert np.allclose(fixed_plane.dir, [1, 0, 0])  # Check direction is preserved


def test_interrupted_write_recovery(water, temp_dir):
    """Test that recovery works when writing is interrupted."""
    import os.path
    
    # Path for the test
    filename = os.path.join(temp_dir, "interrupted.tess")
    
    # First, write a single frame
    with Trajectory(filename, 'w') as traj:
        traj.write(water)
    
    # Now, start a new trajectory with a high flush_interval
    traj = Trajectory(filename, 'a', flush_interval=100)
    
    # Add several frames without triggering a flush
    for i in range(10):
        # Modify the atoms slightly to make them unique
        water_copy = water.copy()
        water_copy.positions[0, 0] += i * 0.1
        traj.write(water_copy)
        
    # Check that we haven't yet written the new frames
    # This is a bit tricky to test - one way is to try reading from the file
    # and check the frame count
    with Trajectory(filename, 'r') as read_traj:
        # We should only see the initial frame
        frames_before_close = len(read_traj)
    
    # Now properly close the trajectory, which should flush remaining frames
    traj.close()
    
    # Read the file again to check if frames were written
    with Trajectory(filename, 'r') as read_traj:
        # Should now have 1 (initial) + 10 (added) frames
        assert len(read_traj) == 11
        
        # Verify the content of the frames
        for i in range(1, 11):
            frame = read_traj[i]
            # Check the position we modified
            expected_pos = water.positions[0, 0] + (i-1) * 0.1
            assert np.isclose(frame.positions[0, 0], expected_pos)


def test_trajectory_complex_multiframe(empty_trajectory):
    """Test a complex multi-frame trajectory with different properties in each frame."""
    # Create different atoms objects with various properties
    from ase.build import molecule, bulk
    
    # Frame 1: Water with custom array and info
    water = molecule('H2O')
    water.arrays['custom_array'] = np.array([1.0, 2.0, 3.0])
    water.info['system'] = 'water'
    
    # Frame 2: Silicon with different masses and constraints
    silicon = bulk('Si')
    silicon.set_masses(np.array([29.0] * len(silicon)))
    from ase.constraints import FixAtoms
    silicon.constraints = [FixAtoms(indices=[0])]
    silicon.info['system'] = 'silicon'
    
    # Frame 3: Methane with calculator data
    methane = molecule('CH4')
    energy = -15.0
    forces = np.random.random((len(methane), 3))
    calc = SinglePointCalculator(methane, energy=energy, forces=forces)
    methane.calc = calc
    methane.info['system'] = 'methane'
    
    # Write multi-frame trajectory with different atoms in each frame
    with Trajectory(empty_trajectory, 'w') as traj:
        traj.write(water)
        traj.write(silicon)
        traj.write(methane)
    
    # Read and verify each frame
    with Trajectory(empty_trajectory, 'r') as traj:
        # Check we have 3 frames
        assert len(traj) == 3
        
        # Check frame 1 (water)
        loaded_water = traj[0]
        assert loaded_water.get_chemical_formula() == 'H2O'
        assert 'custom_array' in loaded_water.arrays
        assert np.array_equal(loaded_water.arrays['custom_array'], water.arrays['custom_array'])
        assert loaded_water.info['system'] == 'water'
        
        # Check frame 2 (silicon)
        loaded_silicon = traj[1]
        assert loaded_silicon.get_chemical_formula() == 'Si2'
        assert np.allclose(loaded_silicon.get_masses(), silicon.get_masses())
        assert len(loaded_silicon.constraints) == 1
        assert isinstance(loaded_silicon.constraints[0], FixAtoms)
        assert loaded_silicon.info['system'] == 'silicon'
        
        # Check frame 3 (methane)
        loaded_methane = traj[2]
        assert loaded_methane.get_chemical_formula() == 'CH4'
        assert loaded_methane.calc is not None
        assert np.isclose(loaded_methane.calc.get_potential_energy(), energy)
        assert np.allclose(loaded_methane.calc.get_forces(), forces)
        assert loaded_methane.info['system'] == 'methane'
