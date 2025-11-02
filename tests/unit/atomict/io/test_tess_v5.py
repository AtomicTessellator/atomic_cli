import tempfile
from pathlib import Path

import numpy as np
import pytest


# Skip if msgpack dependencies are missing
pytest.importorskip("msgpack")
pytest.importorskip("msgpack_numpy")

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from atomict.io.formats.tess import write_tess, read_tess


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tdir:
        yield Path(tdir)


def _make_base_atoms() -> Atoms:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.74],
        ],
        dtype=np.float64,
    )
    cell = np.diag([5.0, 5.0, 5.0]).astype(np.float64)
    atoms = Atoms(symbols="H2", positions=positions, pbc=True)
    atoms.set_cell(cell)
    return atoms


def _set_calc_results(atoms: Atoms, energy: float, forces: np.ndarray, stress: np.ndarray) -> None:
    calc = SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress)
    atoms.calc = calc


def test_v5_single_frame_full_state_roundtrip(temp_dir):
    a = _make_base_atoms()

    # Static arrays
    a.set_tags(np.array([1, 2], dtype=int))
    a.set_initial_charges(np.array([0.1, -0.1], dtype=float))
    a.set_initial_magnetic_moments(np.array([0.0, 1.0], dtype=float))

    # Per-frame arrays
    momenta = np.array([[0.1, 0.0, 0.0], [0.0, -0.1, 0.0]], dtype=float)
    a.set_momenta(momenta)

    # Forces both as array and in calculator
    forces = np.array([[0.01, 0.0, 0.0], [0.0, -0.02, 0.0]], dtype=float)
    a.arrays['forces'] = forces.copy()

    # Stress and energy via calculator
    stress = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03], dtype=float)
    energy = -1.234
    _set_calc_results(a, energy=energy, forces=forces, stress=stress)

    # Info and custom arrays
    a.info['source'] = 'unit-test'
    a.info['int_val'] = 7
    a.arrays['my_custom'] = np.array([10.0, 20.0], dtype=float)

    out_path = temp_dir / "full_v5_single.tess"
    write_tess(a, str(out_path), metadata={"k": "v"}, compression="zlib", compression_level=1)

    out_atoms, out_meta = read_tess(str(out_path))
    assert len(out_atoms) == 1
    b = out_atoms[0]

    # Core state
    assert a.get_chemical_symbols() == b.get_chemical_symbols()
    assert np.allclose(a.get_positions(), b.get_positions())
    assert np.allclose(a.get_cell(), b.get_cell())
    assert (a.get_pbc() == b.get_pbc()).all()
    assert np.allclose(a.get_masses(), b.get_masses())

    # Static arrays and per-frame arrays
    assert np.allclose(a.get_tags(), b.get_tags())
    assert np.allclose(a.get_initial_charges(), b.get_initial_charges())
    assert np.allclose(a.get_initial_magnetic_moments(), b.get_initial_magnetic_moments())
    assert np.allclose(a.get_momenta(), b.get_momenta())

    # Forces array
    assert 'forces' in b.arrays
    assert np.allclose(forces, b.arrays['forces'])

    # Calculator was reconstructed with results
    assert b.calc is not None
    assert b.calc.results
    assert pytest.approx(energy) == b.calc.results.get('energy')
    assert np.allclose(forces, b.calc.results.get('forces'))
    assert np.allclose(stress, b.calc.results.get('stress'))

    # Info and custom arrays
    assert b.info.get('source') == 'unit-test'
    assert b.info.get('int_val') == 7
    assert 'my_custom' in b.arrays and np.allclose(b.arrays['my_custom'], np.array([10.0, 20.0]))
    assert out_meta == {"k": "v"}


def test_v5_multi_frame_static_vs_per_frame_arrays(temp_dir):
    frames = []
    base = _make_base_atoms()
    base.set_tags(np.array([1, 1], dtype=int))
    base.set_initial_charges(np.array([0.2, -0.2], dtype=float))
    base.set_initial_magnetic_moments(np.array([0.5, 0.5], dtype=float))

    for i in range(4):
        f = base.copy()
        # vary positions deterministically
        f.set_positions(f.get_positions() + np.array([0.0, 0.0, 0.05 * i]))
        # vary momenta per-frame
        f.set_momenta(np.array([[i*0.1, 0.0, 0.0], [0.0, -i*0.1, 0.0]], dtype=float))
        # per-frame forces present only on even frames
        if i % 2 == 0:
            F = np.array([[0.01*i, 0.0, 0.0], [0.0, -0.02*i, 0.0]], dtype=float)
            f.arrays['forces'] = F
        # per-frame stress via calc (simulate)
        stress = np.array([0.1*i, 0.2*i, 0.3*i, 0.01*i, 0.02*i, 0.03*i], dtype=float)
        E = -1.0 - i*0.1
        _set_calc_results(f, energy=E, forces=f.arrays.get('forces'), stress=stress)
        frames.append(f)

    out_path = temp_dir / "v5_static_perframe.tess"
    write_tess(frames, str(out_path), compression="zlib", compression_level=1)

    out_atoms, _ = read_tess(str(out_path))
    assert len(out_atoms) == len(frames)

    for i, (a, b) in enumerate(zip(frames, out_atoms)):
        # invariants
        assert np.allclose(a.get_tags(), b.get_tags())
        assert np.allclose(a.get_initial_charges(), b.get_initial_charges())
        assert np.allclose(a.get_initial_magnetic_moments(), b.get_initial_magnetic_moments())
        # per-frame
        assert np.allclose(a.get_positions(), b.get_positions())
        assert np.allclose(a.get_momenta(), b.get_momenta())
        # forces: may be missing on odd frames
        if 'forces' in a.arrays:
            assert 'forces' in b.arrays
            assert np.allclose(a.arrays['forces'], b.arrays['forces'])
        else:
            assert 'forces' not in b.arrays or b.arrays['forces'] is None or b.arrays['forces'].size == 0
        # calculator present and contains stress/energy for each frame
        assert b.calc is not None and b.calc.results
        assert 'stress' in b.calc.results
        assert 'energy' in b.calc.results


def test_v5_subset_frames_preserves_properties(temp_dir):
    frames = []
    for i in range(5):
        a = _make_base_atoms()
        a.set_momenta(np.array([[i, 0.0, 0.0], [0.0, -i, 0.0]], dtype=float))
        a.arrays['forces'] = np.array([[i*0.01, 0.0, 0.0], [0.0, -i*0.02, 0.0]], dtype=float)
        stress = np.array([i*0.1, 0, 0, 0, 0, 0], dtype=float)
        _set_calc_results(a, energy=-i, forces=a.arrays['forces'], stress=stress)
        a.info['i'] = i
        frames.append(a)

    out_path = temp_dir / "v5_subset.tess"
    write_tess(frames, str(out_path), compression="zlib", compression_level=1)

    subset = [4, 1, 3]
    out_atoms, _ = read_tess(str(out_path), frames_indices=subset)
    assert [at.info.get('i') for at in out_atoms] == subset
    for idx, at in zip(subset, out_atoms):
        assert np.allclose(at.get_momenta()[0, 0], float(idx))
        assert np.allclose(at.calc.results['energy'], -float(idx))


def test_v5_variable_species_raises(temp_dir):
    a = _make_base_atoms()
    b = Atoms(symbols="H3", positions=np.vstack([a.get_positions(), [[0.1, 0.1, 0.1]]]), pbc=True)
    b.set_cell(a.get_cell())
    with pytest.raises(ValueError):
        write_tess([a, b], str(temp_dir / "invalid_species.tess"))


