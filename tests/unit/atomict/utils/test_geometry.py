import pytest
import numpy as np
from atomict.utils.fhiaims.geometry import fhi_to_ase
from ase import Atoms


@pytest.fixture
def periodic_geometry_file(tmp_path):
    content = """
lattice_vector 5.0 0.0 0.0
lattice_vector 0.0 5.0 0.0
lattice_vector 0.0 0.0 5.0
atom 0.0 0.0 0.0 Si
atom 2.5 2.5 2.5 Si
    """.strip()
    
    file_path = tmp_path / "geometry.in"
    with open(file_path, "w") as f:
        f.write(content)
    return str(file_path)


@pytest.fixture
def non_periodic_geometry_file(tmp_path):
    content = """
atom 0.0 0.0 0.0 H
atom 0.0 0.0 0.74 H
    """.strip()
    
    file_path = tmp_path / "geometry.in"
    with open(file_path, "w") as f:
        f.write(content)
    return str(file_path)


def test_fhi_to_ase_periodic(periodic_geometry_file):
    """Test conversion of periodic system (with lattice vectors)"""
    atoms = fhi_to_ase(periodic_geometry_file)
    
    # Check if it's an ASE Atoms object
    assert isinstance(atoms, Atoms)
    
    # Check periodicity
    assert all(atoms.pbc)
    
    # Check lattice vectors
    expected_cell = np.array([[5.0, 0.0, 0.0],
                            [0.0, 5.0, 0.0],
                            [0.0, 0.0, 5.0]])
    assert np.allclose(atoms.cell, expected_cell)
    
    # Check atoms
    assert len(atoms) == 2
    assert all(symbol == 'Si' for symbol in atoms.symbols)
    assert np.allclose(atoms.positions[0], [0.0, 0.0, 0.0])
    assert np.allclose(atoms.positions[1], [2.5, 2.5, 2.5])


def test_fhi_to_ase_non_periodic(non_periodic_geometry_file):
    """Test conversion of non-periodic system (no lattice vectors)"""
    atoms = fhi_to_ase(non_periodic_geometry_file)
    
    # Check if it's an ASE Atoms object
    assert isinstance(atoms, Atoms)
    
    # Check non-periodicity
    assert not any(atoms.pbc)
    
    # Check atoms
    assert len(atoms) == 2
    assert all(symbol == 'H' for symbol in atoms.symbols)
    assert np.allclose(atoms.positions[0], [0.0, 0.0, 0.0])
    assert np.allclose(atoms.positions[1], [0.0, 0.0, 0.74])


def test_fhi_to_ase_file_not_found():
    """Test handling of non-existent file"""
    with pytest.raises(FileNotFoundError):
        fhi_to_ase("nonexistent_file.in")