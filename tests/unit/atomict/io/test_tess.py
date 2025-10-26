import os
from pathlib import Path
import tempfile
import importlib.util
import pytest


# Skip all tests in this module if msgpack dependencies are missing
pytest.importorskip("msgpack")
pytest.importorskip("msgpack_numpy")

from atomict.io.formats.tess import write_tess, read_tess


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tdir:
        yield Path(tdir)


def _have_lz4() -> bool:
    return importlib.util.find_spec("lz4.block") is not None


def _make_simple_atoms():
    # Create a tiny diatomic system with a simple cubic cell
    from ase import Atoms
    import numpy as np

    positions = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.74],
    ])
    cell = np.diag([5.0, 5.0, 5.0])
    atoms = Atoms(symbols="H2", positions=positions, pbc=True)
    atoms.set_cell(cell)
    return atoms


def _compare_atoms(a, b):
    import numpy as np

    assert a.get_chemical_symbols() == b.get_chemical_symbols()
    assert np.allclose(a.get_positions(), b.get_positions(), rtol=1e-6, atol=1e-6)
    assert np.allclose(a.get_cell(), b.get_cell(), rtol=1e-6, atol=1e-6)
    assert (a.get_pbc() == b.get_pbc()).all()
    # Masses are set from header in v4
    assert np.allclose(a.get_masses(), b.get_masses(), rtol=0, atol=0)


@pytest.mark.parametrize(
    "compression",
    [
        "none",
        "zlib",
        pytest.param("lz4", marks=pytest.mark.skipif(not _have_lz4(), reason="lz4 not installed")),
    ],
)
def test_roundtrip_single_frame(temp_dir, compression):
    atoms = _make_simple_atoms()
    tess_path = temp_dir / "single.tess"

    meta = {"creator": "unit-test", "version": 1}
    write_tess(atoms, str(tess_path), metadata=meta, compression=compression, compression_level=1)

    out_atoms, out_meta = read_tess(str(tess_path))
    assert isinstance(out_atoms, list) and len(out_atoms) == 1
    _compare_atoms(atoms, out_atoms[0])
    assert out_meta == meta


def test_roundtrip_multi_frame_same_cell(temp_dir):
    from ase import Atoms
    import numpy as np

    base = _make_simple_atoms()
    frames = []
    for i in range(5):
        f = base.copy()
        shift = np.array([0.0, 0.0, 0.05 * i])
        f.set_positions(f.get_positions() + shift)
        frames.append(f)

    tess_path = temp_dir / "same_cell.tess"
    write_tess(frames, str(tess_path), compression="zlib", compression_level=1)

    out_atoms, _ = read_tess(str(tess_path))
    assert len(out_atoms) == len(frames)
    for a, b in zip(frames, out_atoms):
        _compare_atoms(a, b)


def test_roundtrip_multi_frame_changing_cell(temp_dir):
    import numpy as np

    frames = []
    for i in range(5):
        a = _make_simple_atoms()
        cell = np.diag([5.0 + 0.1 * i, 5.0, 5.0])
        a.set_cell(cell)
        frames.append(a)

    tess_path = temp_dir / "changing_cell.tess"
    write_tess(frames, str(tess_path), compression="zlib", compression_level=1)

    out_atoms, _ = read_tess(str(tess_path))
    assert len(out_atoms) == len(frames)
    for a, b in zip(frames, out_atoms):
        _compare_atoms(a, b)


def test_read_subset_frames_with_indices(temp_dir):
    import numpy as np

    frames = []
    for i in range(6):
        a = _make_simple_atoms()
        # encode index in z position for easy checking
        pos = a.get_positions()
        pos[0, 2] = float(i)
        a.set_positions(pos)
        frames.append(a)

    tess_path = temp_dir / "subset.tess"
    write_tess(frames, str(tess_path), compression="zlib", compression_level=1)

    requested = [4, 1, 5]
    out_atoms, _ = read_tess(str(tess_path), frames_indices=requested)
    assert len(out_atoms) == len(requested)
    for expected_idx, atoms in zip(requested, out_atoms):
        assert atoms.get_positions()[0, 2] == pytest.approx(float(expected_idx))


def test_invalid_compression_raises(temp_dir):
    atoms = _make_simple_atoms()
    tess_path = temp_dir / "invalid.tess"
    with pytest.raises(ValueError):
        write_tess(atoms, str(tess_path), compression="gzip")


def test_read_out_of_range_index_raises(temp_dir):
    frames = [_make_simple_atoms() for _ in range(2)]
    tess_path = temp_dir / "oor.tess"
    write_tess(frames, str(tess_path), compression="none")

    with pytest.raises(IndexError):
        read_tess(str(tess_path), frames_indices=[2])

