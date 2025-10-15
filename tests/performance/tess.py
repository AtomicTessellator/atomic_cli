import os
import sys
import math
import time
from typing import List, Tuple


def _resolve_repo_root() -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    # tests/performance -> tests -> repo_root
    return os.path.abspath(os.path.join(here, '..', '..'))


def _load_helpers():
    repo_root = _resolve_repo_root()
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    # Reuse helpers from format_test.py
    import format_test as ft  # type: ignore
    return ft


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
    return supercell, (a, b, c)


def _make_frames_with_jitter(template, num_frames: int, jitter_sigma: float = 0.02, seed: int = 0):
    import numpy as np

    rng = np.random.default_rng(seed)
    frames = []
    for _ in range(num_frames):
        fr = template.copy()
        pos = fr.get_positions()
        pos += rng.normal(0.0, jitter_sigma, size=pos.shape)
        fr.set_positions(pos)
        frames.append(fr)
    return frames


def _write_traj(frames, out_path: str) -> float:
    from atomict.io.formats.traj import write_traj

    start = time.perf_counter()
    write_traj(frames, out_path, metadata=None)
    return time.perf_counter() - start


def _read_traj(in_path: str) -> Tuple[List['ase.Atoms'], dict, float]:
    from atomict.io.formats.traj import read_traj

    start = time.perf_counter()
    atoms_list, metadata = read_traj(in_path)
    return atoms_list, metadata, time.perf_counter() - start


def _write_tess(frames, out_path: str) -> float:
    from atomict.io.msgpack import save_msgpack_trajectory

    start = time.perf_counter()
    save_msgpack_trajectory(frames, out_path, metadata=None)
    return time.perf_counter() - start


def _read_tess(in_path: str) -> Tuple[List['ase.Atoms'], dict, float]:
    from atomict.io.msgpack import load_msgpack_trajectory

    start = time.perf_counter()
    atoms_list, metadata = load_msgpack_trajectory(in_path)
    return atoms_list, metadata, time.perf_counter() - start


def main():
    ft = _load_helpers()

    repo_root = _resolve_repo_root()
    cif_path = os.path.join(repo_root, 'tests', 'fixtures', 'Al.cif')
    if not os.path.exists(cif_path):
        raise FileNotFoundError(cif_path)

    # Load base structure
    from ase.io import read as ase_read

    unit = ase_read(cif_path)

    # Build ~100_000-atom supercell
    target_atoms = 100_000
    supercell, repeat = _build_supercell_to_target_atoms(unit, target_atoms)
    natoms = len(supercell)

    # Create 100 jittered frames
    nframes = 100
    frames = _make_frames_with_jitter(supercell, nframes, jitter_sigma=0.02, seed=0)

    # Output paths in /tmp
    base = f"/tmp/atomict_perf.Al_{natoms}x{nframes}"
    traj_path = base + '.traj'
    tess_path = base + '.tess'

    # Ensure parent exists
    os.makedirs(os.path.dirname(traj_path), exist_ok=True)

    # Write
    traj_write_s = _write_traj(frames, traj_path)
    tess_write_s = _write_tess(frames, tess_path)

    # Sizes
    traj_mb = ft.get_file_size_mb(traj_path)
    tess_mb = ft.get_file_size_mb(tess_path)

    # Read back
    traj_atoms, traj_meta, traj_read_s = _read_traj(traj_path)
    tess_atoms, tess_meta, tess_read_s = _read_tess(tess_path)

    # Report
    print(f"Supercell repeats: {repeat} -> atoms/frame={natoms}, frames={nframes}")
    header = (
        ("FORMAT", 8),
        ("WRITE(s)", 10),
        ("READ(s)", 10),
        ("SIZE(MB)", 10),
        ("FRAMES", 8),
        ("ATOMS/FR", 9),
    )
    line = ''.join(h[:w].ljust(w) for h, w in header)
    sep = '-' * len(line)
    print(sep)
    print(line)
    print(sep)

    def row(fmt: str, w_s: float, r_s: float, size_mb: float):
        print(
            f"{fmt:<{header[0][1]}.{header[0][1]}}"
            f"{w_s:>{header[1][1]}.3f}"
            f"{r_s:>{header[2][1]}.3f}"
            f"{size_mb:>{header[3][1]}.2f}"
            f"{nframes:>{header[4][1]}}"
            f"{natoms:>{header[5][1]}}"
        )

    row('traj', traj_write_s, traj_read_s, traj_mb)
    row('tess', tess_write_s, tess_read_s, tess_mb)
    print(sep)


if __name__ == '__main__':
    main()


