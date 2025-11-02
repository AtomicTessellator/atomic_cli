import math
from pathlib import Path

import numpy as np
import pytest


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _build_supercell_to_target_atoms(unit_cell, target_atoms: int):
    n_base = len(unit_cell)
    if n_base <= 0:
        raise ValueError("CIF parsed 0 atoms")
    mult_total = int(math.ceil(target_atoms / float(n_base)))
    side = int(math.ceil(mult_total ** (1.0 / 3.0)))
    a = side
    b = side
    c = int(math.ceil(mult_total / float(a * b)))
    supercell = unit_cell.repeat((a, b, c))
    return supercell


def _make_frames(template, num_frames: int, jitter_sigma: float = 0.01, seed: int = 0):
    rng = np.random.default_rng(seed)
    frames = []
    for _ in range(num_frames):
        fr = template.copy()
        pos = fr.get_positions()
        pos += rng.normal(0.0, jitter_sigma, size=pos.shape)
        fr.set_positions(pos)
        frames.append(fr)
    return frames


@pytest.mark.slow
def test_large_atoms_per_frame_tess(tmp_path):
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from ase.io import read as ase_read
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    cif_path = _repo_root() / 'tests' / 'fixtures' / 'Al.cif'
    unit = ase_read(str(cif_path))

    supercell = _build_supercell_to_target_atoms(unit, target_atoms=3000)
    frames = _make_frames(supercell, num_frames=3)

    out_path = tmp_path / 'large.tess'
    with TrajectoryWriter(str(out_path), mode='w', tess_compression='zlib', tess_compression_level=1) as tw:
        for fr in frames:
            tw.write(fr)

    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) == len(frames)
        # Validate atom counts
        for fr in tr:
            assert len(fr) == len(supercell)


@pytest.mark.slow
def test_large_atoms_per_frame_atraj(tmp_path):
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    from ase.io import read as ase_read
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    cif_path = _repo_root() / 'tests' / 'fixtures' / 'Al.cif'
    unit = ase_read(str(cif_path))

    supercell = _build_supercell_to_target_atoms(unit, target_atoms=2500)
    frames = _make_frames(supercell, num_frames=2)

    out_path = tmp_path / 'large.atraj'
    with TrajectoryWriter(str(out_path), mode='w') as tw:
        for fr in frames:
            tw.write(fr)

    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) == len(frames)
        for fr in tr:
            assert len(fr) == len(supercell)


