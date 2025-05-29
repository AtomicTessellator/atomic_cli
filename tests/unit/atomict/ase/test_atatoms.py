import numpy as np
from ase.build import molecule, bulk
from ase import Atom
from ase.calculators.emt import EMT
from ase.optimize import BFGS
from atomict.ase.atatoms import ATAtoms
from ase.atoms import Atoms


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


##########################################
### below tests are for manual testing ###
##########################################
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


def create_malachite():
    # Malachite crystallizes in the monoclinic system
    # Space group: P21/a (No. 14)
    
    # Lattice parameters (in Angstroms and degrees)
    a = 9.502
    b = 11.974
    c = 3.240
    alpha = 90.0
    beta = 98.75
    gamma = 90.0
    
    # Fractional coordinates of atoms in the unit cell
    # These are approximate positions based on crystallographic data
    # Cu positions
    cu1_pos = (0.0, 0.0, 0.0)
    cu2_pos = (0.5, 0.5, 0.0)
    
    # C position
    c_pos = (0.25, 0.25, 0.5)
    
    # O positions (including both carbonate oxygens and hydroxyl oxygens)
    o1_pos = (0.3, 0.3, 0.5)   # carbonate oxygen
    o2_pos = (0.2, 0.15, 0.5)  # carbonate oxygen
    o3_pos = (0.25, 0.35, 0.5) # carbonate oxygen
    o4_pos = (0.1, 0.05, 0.0)  # hydroxyl oxygen
    o5_pos = (0.4, 0.45, 0.0)  # hydroxyl oxygen
    
    # H positions (for hydroxyl groups)
    h1_pos = (0.15, 0.08, 0.0)
    h2_pos = (0.35, 0.42, 0.0)
    
    # Create the atoms object with the unit cell
    malachite = Atoms('Cu2CO5H2',
                     positions=[cu1_pos, cu2_pos, c_pos, 
                               o1_pos, o2_pos, o3_pos, o4_pos, o5_pos, 
                               h1_pos, h2_pos],
                     cell=[a, b, c, alpha, beta, gamma],
                     pbc=True)
    
    # Convert fractional to cartesian coordinates
    cell = malachite.get_cell()
    positions = malachite.get_positions()
    positions = np.dot(positions, cell)
    malachite.set_positions(positions)
    
    return malachite


def test_context_manager_cleanup():
    """Test that the context manager properly cleans up resources."""
    import threading
    from unittest.mock import Mock, patch
    import time
    import queue
    
    # Create a mock worker that simulates the real worker's behavior
    mock_worker = Mock()
    mock_worker._queue = queue.Queue()
    
    # Create a thread that will stay alive until explicitly stopped
    stop_event = threading.Event()
    def worker_thread():
        while not stop_event.is_set():
            try:
                # Simulate processing the queue
                item = mock_worker._queue.get(timeout=0.1)
                if item is None:  # Signal to stop
                    break
                mock_worker._queue.task_done()
            except queue.Empty:
                pass
    
    mock_worker._thread = threading.Thread(target=worker_thread)
    mock_worker._thread.start()

    try:
        # Create a mock ATAtoms instance
        with patch('atomict.ase.atatoms.ATAtoms._create_worker', return_value=mock_worker) as mock_create_worker:
            with ATAtoms(create_malachite(), batch_diffs=True) as at_atoms:
                # Verify worker was created
                mock_create_worker.assert_called_once()
                
                # Verify thread is alive during context
                assert mock_worker._thread.is_alive()
                
                # Simulate some operations
                original_positions = at_atoms.get_positions().copy()
                at_atoms.set_positions(original_positions)
                at_atoms.set_positions(original_positions * 0.8)  # 20% compression
        
        # Verify thread is no longer alive after context
        assert not mock_worker._thread.is_alive()
    finally:
        # Ensure the thread is stopped even if the test fails
        stop_event.set()
        mock_worker._thread.join(timeout=1.0)
