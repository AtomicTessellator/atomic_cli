import os
import numpy as np
from ase import Atom
from ase.calculators.emt import EMT
from ase.calculators.aims import Aims
from ase.optimize import BFGS
from atomict.ase.atatoms import ATAtoms
from ase.atoms import Atoms
from ase.filters import UnitCellFilter
import pytest
from unittest.mock import patch, Mock

from tests.unit.atomict.ase.test_atatoms import create_malachite


@pytest.fixture
def mock_post():
    """Mock the post function from atomict.api as imported in atatoms module"""
    with patch('atomict.ase.atatoms.post') as mock:
        mock.return_value = {'id': 'test-state-id'}
        yield mock


@pytest.fixture
def mock_patch():
    """Mock the patch function from atomict.api as imported in atatoms module"""
    with patch('atomict.ase.atatoms.patch') as mock:
        mock.return_value = {'id': 'test-run-id'}
        yield mock


@pytest.fixture
def mock_api_calls(mock_post, mock_patch):
    """Fixture that combines both API mocks"""
    return {'post': mock_post, 'patch': mock_patch}


def test_atatoms_integration(mock_api_calls):
    """Test all ASE atoms operations on ATAtoms objects."""
    # Initialize ATAtoms with project_id
    malachite = create_malachite()
    project_id = os.environ.get('AT_PROJECT_ID', 'test-project-id')
    atoms = ATAtoms(malachite, project_id=project_id)
    
    # Geometry optimization
    atoms.calc = EMT()

    # aims calc
    # from ase.calculators.aims import AimsProfile
    # import os
    # species_dir = '/home/steven/test/species_defaults/defaults_2020/light'
    # profile = AimsProfile(os.environ['AIMS_COMMAND'], default_species_dir=species_dir)
    # calc = Aims(
    #     profile=profile,
    #     xc='pbe0',
    #     species_dir=species_dir,
    #     k_grid=(4, 4, 4)  # Add k-point sampling
    # )
    # atoms.calc = calc
    # opt = BFGS(atoms)
    # opt.run(fmax=0.01)

    # forces
    original_positions = atoms.get_positions().copy()
    atoms.set_positions(original_positions)
    atoms.set_positions(original_positions * 0.8)  # 20% compression

    # Cell relaxation (shape/volume change)
    atoms.calc = EMT()
    ucf = UnitCellFilter(atoms)
    opt = BFGS(ucf)
    opt.run(fmax=0.02)
    
    # Isotropic rescaling - expands atoms
    original_cell = atoms.get_cell().array.copy()
    atoms.set_cell(atoms.get_cell() * 4.00, scale_atoms=True)
    assert np.allclose(atoms.get_cell().array, original_cell * 4.00)
    
    # Affine strain application - stretch cell along x, y, z axis by %
    current_cell = atoms.get_cell().array.copy()
    strain = np.array([[1.50, 0, 0], [0, 1.25, 0], [0, 0, 0.66]])
    atoms.set_cell(np.dot(atoms.get_cell(), strain), scale_atoms=True)
    assert not np.allclose(atoms.get_cell().array, current_cell)
    
    # Adding vacuum
    pre_vacuum_volume = atoms.get_cell().volume
    atoms.center(vacuum=10.0)
    assert atoms.get_cell().volume > pre_vacuum_volume
    
    # Manual atomic displacement
    original_pos = atoms[0].position.copy()
    atoms[0].position += [10, 10.0, 10.0]
    assert np.allclose(atoms[0].position, original_pos + [10, 10.0, 10.0])
    
    # Deleting/adding atoms
    n_atoms = len(atoms)
    del atoms[0]
    assert len(atoms) == n_atoms - 1
    
    atoms += Atoms('Og', positions=[[5.0, 5.0, 5.0]])
    assert len(atoms) == n_atoms
    assert atoms.get_chemical_symbols()[-1] == 'Og'
    
    # Reordering atoms
    orig_symbols = atoms.get_chemical_symbols()
    atoms = atoms[::-1]
    assert atoms.get_chemical_symbols() == orig_symbols[::-1]
    assert isinstance(atoms, ATAtoms)
    
    # Changing atomic charges
    charges = [0.1] * len(atoms)
    atoms.set_initial_charges(charges)
    assert np.allclose(atoms.get_initial_charges(), charges)
    
    # Assigning magnetic moments
    mag_moments = [2.0] * len(atoms)
    atoms.set_initial_magnetic_moments(mag_moments)
    assert np.allclose(atoms.get_initial_magnetic_moments(), mag_moments)
    
    # Changing atomic masses
    masses = [10.0] * len(atoms)
    atoms.set_masses(masses)
    assert np.allclose(atoms.get_masses(), masses)
    
    # Assigning forces manually (e.g., post-processing)
    random_forces = np.random.randn(len(atoms), 3)
    atoms.set_array('forces', random_forces)
    assert np.allclose(atoms.get_array('forces'), random_forces)
    
    # Supercell creation
    n_atoms = len(atoms)
    atoms = atoms.repeat((2, 2, 2))
    assert len(atoms) == n_atoms * 8
    assert isinstance(atoms, ATAtoms)
    
    # Rotation
    pre_rotation_pos = atoms.get_positions().copy()
    atoms.rotate('z', 90, center='COM')
    assert not np.allclose(atoms.get_positions(), pre_rotation_pos)
    
    # Translation
    pre_translation_pos = atoms.get_positions().copy()
    atoms.translate([1.0, 0.0, 0.0])
    assert np.allclose(atoms.get_positions(), pre_translation_pos + [1.0, 0.0, 0.0])
    
    # Random thermal displacement
    pre_displacement_pos = atoms.get_positions().copy()
    displacements = np.random.normal(0, 0.05, size=(len(atoms), 3))
    atoms.set_positions(atoms.get_positions() + displacements)
    assert not np.allclose(atoms.get_positions(), pre_displacement_pos)
    
    # Save final state
    state_id = atoms.save_current_state()
    assert state_id is not None