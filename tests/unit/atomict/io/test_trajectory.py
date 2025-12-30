import tempfile
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tdir:
        yield Path(tdir)


def _make_atoms() -> 'ase.Atoms':
    from ase import Atoms

    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.74],
    ])
    cell = np.diag([5.0, 5.0, 5.0])
    atoms = Atoms(symbols="H2", positions=positions, pbc=True)
    atoms.set_cell(cell)
    return atoms


def _compare_atoms(a, b):
    assert a.get_chemical_symbols() == b.get_chemical_symbols()
    assert np.allclose(a.get_positions(), b.get_positions(), rtol=1e-6, atol=1e-6)
    assert np.allclose(a.get_cell(), b.get_cell(), rtol=1e-6, atol=1e-6)
    assert (a.get_pbc() == b.get_pbc()).all()


def test_roundtrip_traj_backend(temp_dir):
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    frames = []
    base = _make_atoms()
    for i in range(3):
        f = base.copy()
        f.set_positions(f.get_positions() + np.array([0.0, 0.0, 0.05 * i]))
        frames.append(f)

    out_path = temp_dir / "ase.traj"
    meta = {"creator": "unit-test", "k": 1}

    with TrajectoryWriter(str(out_path), mode='w', metadata=meta) as tw:
        for fr in frames:
            tw.write(fr)

    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) == len(frames)
        for a, b in zip(frames, tr):
            _compare_atoms(a, b)


def test_roundtrip_atraj_backend(temp_dir):
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    frames = []
    base = _make_atoms()
    for i in range(4):
        f = base.copy()
        f.set_positions(f.get_positions() + np.array([0.0, 0.0, 0.03 * i]))
        frames.append(f)

    out_path = temp_dir / "test.atraj"
    meta = {"project": "traj-unit"}

    with TrajectoryWriter(str(out_path), mode='w', metadata=meta) as tw:
        for fr in frames:
            tw.write(fr)

    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) == len(frames)
        # metadata property should be available
        assert isinstance(tr.metadata, dict)
        for a, b in zip(frames, tr):
            _compare_atoms(a, b)


def test_roundtrip_tess_backend(temp_dir):
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    frames = []
    base = _make_atoms()
    for i in range(5):
        f = base.copy()
        f.set_positions(f.get_positions() + np.array([0.0, 0.0, 0.02 * i]))
        frames.append(f)

    out_path = temp_dir / "test.tess"
    meta = {"source": "unittest"}

    with TrajectoryWriter(str(out_path), mode='w', metadata=meta, tess_compression="zlib", tess_compression_level=1) as tw:
        for fr in frames:
            tw.write(fr)

    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) == len(frames)
        for idx in range(len(frames)):
            _compare_atoms(frames[idx], tr[idx])
        # iterate
        assert sum(1 for _ in tr) == len(frames)


def test_calculator_properties_preserved(temp_dir):
    """Test that calculator results (energy, forces, stress) are preserved when buffering atoms."""
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from ase import Atoms
    from ase.calculators.singlepoint import SinglePointCalculator
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    # Create atoms with calculator results
    frames_with_calc = []
    base = _make_atoms()
    
    for i in range(3):
        atoms = base.copy()
        atoms.set_positions(atoms.get_positions() + np.array([0.0, 0.0, 0.05 * i]))
        
        # Attach a SinglePointCalculator with mock results
        energy = -10.0 - i * 0.5
        forces = np.array([
            [0.1 * i, 0.0, -0.2],
            [-0.1 * i, 0.0, 0.2],
        ])
        stress = np.array([0.01, 0.02, 0.03, 0.0, 0.0, 0.0]) * (i + 1)
        
        calc = SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress)
        atoms.calc = calc
        
        frames_with_calc.append(atoms)
    
    # Test with .atraj format (uses buffering, not ASE backend)
    out_path = temp_dir / "calc_test.atraj"
    
    with TrajectoryWriter(str(out_path), mode='w') as tw:
        for fr in frames_with_calc:
            tw.write(fr)
    
    # Read back and verify calculator properties are accessible
    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) == len(frames_with_calc)
        
        for i, (original, loaded) in enumerate(zip(frames_with_calc, tr)):
            # Verify atoms match
            _compare_atoms(original, loaded)
            
            # Verify calculator results are preserved
            assert loaded.calc is not None, f"Frame {i}: calculator should be attached"
            
            # Check energy
            loaded_energy = loaded.get_potential_energy()
            original_energy = original.get_potential_energy()
            assert np.isclose(loaded_energy, original_energy, rtol=1e-10), \
                f"Frame {i}: energy mismatch {loaded_energy} vs {original_energy}"
            
            # Check forces
            loaded_forces = loaded.get_forces()
            original_forces = original.get_forces()
            assert np.allclose(loaded_forces, original_forces, rtol=1e-10), \
                f"Frame {i}: forces mismatch"
            
            # Check stress
            loaded_stress = loaded.get_stress()
            original_stress = original.get_stress()
            assert np.allclose(loaded_stress, original_stress, rtol=1e-10), \
                f"Frame {i}: stress mismatch"


def test_calculator_properties_preserved_tess(temp_dir):
    """Test that calculator results are preserved in .tess format."""
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from ase.calculators.singlepoint import SinglePointCalculator
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    # Create atoms with calculator results
    atoms = _make_atoms()
    energy = -5.25
    forces = np.array([
        [0.1, -0.2, 0.3],
        [-0.1, 0.2, -0.3],
    ])
    
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
    atoms.calc = calc
    
    out_path = temp_dir / "calc_test.tess"
    
    with TrajectoryWriter(str(out_path), mode='w') as tw:
        tw.write(atoms)
    
    with TrajectoryReader(str(out_path)) as tr:
        loaded = tr[0]
        
        _compare_atoms(atoms, loaded)
        
        assert loaded.calc is not None, "Calculator should be attached"
        assert np.isclose(loaded.get_potential_energy(), energy, rtol=1e-10)
        assert np.allclose(loaded.get_forces(), forces, rtol=1e-10)


def test_atoms_without_calculator(temp_dir):
    """Test that atoms without calculators still work correctly."""
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    # Create atoms without any calculator
    atoms = _make_atoms()
    assert atoms.calc is None
    
    out_path = temp_dir / "no_calc.atraj"
    
    with TrajectoryWriter(str(out_path), mode='w') as tw:
        tw.write(atoms)
    
    with TrajectoryReader(str(out_path)) as tr:
        loaded = tr[0]
        _compare_atoms(atoms, loaded)
        # No calculator should be attached
        assert loaded.calc is None

