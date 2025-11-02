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


