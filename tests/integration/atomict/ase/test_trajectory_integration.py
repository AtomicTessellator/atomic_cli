import tempfile
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def temp_dir():
    with tempfile.TemporaryDirectory() as tdir:
        yield Path(tdir)


def _make_atoms() -> 'ase.Atoms':
    from ase.build import molecule
    a = molecule('H2')
    a.set_cell(np.diag([5.0, 5.0, 5.0]))
    a.center()
    return a


def test_bfgs_writes_trajectory_traj(temp_dir):
    pytest.importorskip("ase")
    from ase.optimize import BFGS
    from ase.calculators.emt import EMT
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    atoms = _make_atoms()
    atoms.calc = EMT()

    out_path = temp_dir / 'opt.traj'
    with TrajectoryWriter(str(out_path), mode='w') as tw:
        opt = BFGS(atoms, trajectory=tw, logfile=None)
        opt.run(fmax=0.1, steps=5)

    with TrajectoryReader(str(out_path)) as tr:
        # Should contain at least initial + some steps
        assert len(tr) >= 1


def test_bfgs_writes_trajectory_tess(temp_dir):
    pytest.importorskip("msgpack")
    pytest.importorskip("msgpack_numpy")
    pytest.importorskip("ase")
    from ase.optimize import BFGS
    from ase.calculators.emt import EMT
    from atomict.io.trajectory import TrajectoryWriter, TrajectoryReader

    atoms = _make_atoms()
    atoms.calc = EMT()

    out_path = temp_dir / 'opt.tess'
    with TrajectoryWriter(str(out_path), mode='w', tess_compression='zlib', tess_compression_level=1) as tw:
        opt = BFGS(atoms, trajectory=tw, logfile=None)
        opt.run(fmax=0.1, steps=5)

    with TrajectoryReader(str(out_path)) as tr:
        assert len(tr) >= 1


