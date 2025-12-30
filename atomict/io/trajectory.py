from typing import Any, Dict, Iterator, List, Optional, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

try:
    from ase.parallel import world
except Exception:
    # Fallback if ASE parallel is unavailable
    class _DummyWorld:
        rank = 0

        def sum_scalar(self, x):
            return x

    world = _DummyWorld()  # type: ignore[assignment]


class TrajectoryWriter:
    """Unified trajectory writer supporting .traj (ASE), .atraj and .tess formats.

    Mirrors the ASE TrajectoryWriter API for use in optimizers (e.g., BFGS(trajectory=...)).
    Frames are written incrementally for .traj and buffered for .atraj/.tess (flushed on close).
    """

    def __init__(
        self,
        filename: str,
        mode: str = "w",
        atoms: Optional['Atoms'] = None,
        properties: Optional[List[str]] = None,
        master: Optional[bool] = None,
        comm: Any = world,
        metadata: Optional[Dict] = None,
        tess_compression: Optional[str] = "zlib",
        tess_compression_level: int = 1,
    ) -> None:
        """Create a trajectory writer.

        Parameters mirror ASE's writer; extra kwargs configure tess compression and metadata.
        """
        # Deferred imports to avoid hard dependency during import time
        from ase import Atoms  # type: ignore

        if mode not in {"w", "a"}:
            raise ValueError('mode must be "w" or "a".')

        self.filename = filename
        self.mode = mode
        self.atoms_default: Optional[Atoms] = atoms
        self.properties = properties
        self.comm = comm
        self.master = (comm.rank == 0) if master is None else bool(master)
        self._description: Dict[str, Any] = dict(metadata or {})

        # Backend selection by extension
        lower = filename.lower()
        if lower.endswith(".traj"):
            # Stream via ASE backend
            from ase.io.trajectory import TrajectoryWriter as _ASEWriter

            self._backend = _ASEWriter(filename, mode=mode, atoms=atoms, properties=properties, master=self.master, comm=comm)
            self._backend_is_ase = True
            self._buffer: Optional[List[Atoms]] = None
        else:
            # Buffer frames and flush on close for msgpack-based formats
            self._backend = None
            self._backend_is_ase = False
            self._buffer = []
            # Tess config
            self._tess_compression = tess_compression
            self._tess_compression_level = int(tess_compression_level)

    def __enter__(self) -> "TrajectoryWriter":
        return self

    def __exit__(self, exc_type, exc_value, tb) -> None:
        self.close()

    def set_description(self, description: Dict[str, Any]) -> None:
        # match ASE method name for compatibility
        self._description.update(description or {})

    def write(self, atoms: Optional['Atoms'] = None, **kwargs: Any) -> None:
        """Write a single frame. If atoms is None, uses atoms passed to __init__.

        For .traj, delegates to ASE. For .atraj/.tess, buffers frames and writes on close.
        """
        from ase import Atoms  # type: ignore

        image: Optional[Atoms] = atoms if atoms is not None else self.atoms_default
        if image is None:
            raise ValueError("No Atoms provided to write() and no default atoms set.")

        if self._backend_is_ase:
            # Pass through calculator properties if provided in kwargs, ASE handles it
            self._backend.write(image, **kwargs)
        else:
            # Buffer plain ASE Atoms; handle OptimizableAtoms or wrappers
            from ase import Atoms  # type: ignore

            plain: Optional[Atoms] = None
            if isinstance(image, Atoms):
                plain = image
            else:
                # Try common wrapper escape hatches
                get_atoms = getattr(image, 'get_atoms', None)
                if callable(get_atoms):
                    maybe_atoms = get_atoms()
                    if isinstance(maybe_atoms, Atoms):
                        plain = maybe_atoms
                if plain is None:
                    inner = getattr(image, 'atoms', None)
                    if isinstance(inner, Atoms):
                        plain = inner
            if plain is None:
                # Best-effort reconstruction
                numbers = getattr(image, 'get_atomic_numbers', lambda: None)()
                positions = getattr(image, 'get_positions', lambda: None)()
                cell = getattr(image, 'get_cell', lambda: None)()
                pbc = getattr(image, 'get_pbc', lambda: None)()
                if numbers is None or positions is None:
                    raise TypeError('Unsupported atoms-like object; expected ASE Atoms or compatible wrapper')
                plain = Atoms(numbers=numbers, positions=positions, cell=cell, pbc=pbc)
            # Make a copy to avoid all frames pointing to the same object after in-place modifications
            copied = plain.copy()

            # Preserve calculator results using SinglePointCalculator
            # (ASE's Atoms.copy() intentionally does not copy the calculator)
            if plain.calc is not None:
                from ase.calculators.singlepoint import SinglePointCalculator

                results = {}
                for prop in ['energy', 'forces', 'stress', 'dipole', 'charges', 'magmom', 'magmoms']:
                    try:
                        val = plain.calc.get_property(prop, plain, allow_calculation=False)
                        if val is not None:
                            results[prop] = val
                    except Exception:
                        pass
                if results:
                    copied.calc = SinglePointCalculator(copied, **results)

            self._buffer.append(copied)

    def close(self) -> None:
        """Flush buffered frames and/or close backend."""
        if self._backend_is_ase:
            self._backend.close()
            return

        # Write buffered frames using the chosen msgpack-based format
        if not self.master:
            # Only master writes similar to ASE semantics
            self._buffer = []
            return

        if not self._buffer:
            # Nothing to write
            return

        lower = self.filename.lower()
        frames = self._buffer

        if lower.endswith('.atraj'):
            from atomict.io.formats.atraj import write_atraj

            write_atraj(frames, self.filename, metadata=self._description or None)
        elif lower.endswith('.tess'):
            from atomict.io.formats.tess import write_tess

            write_tess(
                frames,
                self.filename,
                metadata=self._description or None,
                compression=self._tess_compression,
                compression_level=self._tess_compression_level,
            )
        else:
            # Fallback: use ASE .traj if unknown extension
            from ase.io.trajectory import TrajectoryWriter as _ASEWriter

            with _ASEWriter(self.filename, mode=self.mode, properties=self.properties, master=self.master, comm=self.comm) as tw:
                for fr in frames:
                    tw.write(fr)

        # Clear buffer after writing
        self._buffer = []

    def __len__(self) -> int:
        if self._backend_is_ase:
            return self.comm.sum_scalar(len(self._backend))
        return len(self._buffer or [])


class TrajectoryReader:
    """Unified trajectory reader supporting .traj (ASE), .atraj and .tess formats.

    Provides random access (__getitem__), slicing, iteration, and length like ASE's reader.
    """

    def __init__(self, filename: str) -> None:
        self.filename = filename
        self._metadata: Dict[str, Any] = {}
        lower = filename.lower()
        if lower.endswith('.traj'):
            from ase.io.trajectory import TrajectoryReader as _ASEReader

            self._backend = _ASEReader(filename)
            self._backend_is_ase = True
            self._frames: Optional[List['Atoms']] = None
        elif lower.endswith('.atraj'):
            from atomict.io.formats.atraj import read_atraj

            frames, meta = read_atraj(filename)
            self._frames = frames
            self._metadata = meta or {}
            self._backend = None
            self._backend_is_ase = False
        elif lower.endswith('.tess'):
            from atomict.io.formats.tess import read_tess

            frames, meta = read_tess(filename)
            self._frames = frames
            self._metadata = meta or {}
            self._backend = None
            self._backend_is_ase = False
        else:
            # Fallback to ASE reader
            from ase.io.trajectory import TrajectoryReader as _ASEReader

            self._backend = _ASEReader(filename)
            self._backend_is_ase = True
            self._frames = None

    def __enter__(self) -> "TrajectoryReader":
        return self

    def __exit__(self, exc_type, exc_value, tb) -> None:
        self.close()

    def close(self) -> None:
        if self._backend_is_ase and self._backend is not None:
            self._backend.close()

    def __len__(self) -> int:
        if self._backend_is_ase and self._backend is not None:
            return len(self._backend)
        return len(self._frames or [])

    def __getitem__(self, i: Union[int, slice]) -> Union['Atoms', List['Atoms']]:
        if isinstance(i, slice):
            if self._backend_is_ase and self._backend is not None:
                return list(self._backend[i])
            return (self._frames or [])[i]
        # int index
        if self._backend_is_ase and self._backend is not None:
            return self._backend[i]
        return (self._frames or [])[i]

    def __iter__(self) -> Iterator['Atoms']:
        if self._backend_is_ase and self._backend is not None:
            for i in range(len(self._backend)):
                yield self._backend[i]
        else:
            for fr in (self._frames or []):
                yield fr

    @property
    def metadata(self) -> Dict[str, Any]:
        return dict(self._metadata)


def Trajectory(filename: str, mode: str = 'r', atoms: Optional['Atoms'] = None, properties: Optional[List[str]] = None, master: Optional[bool] = None, comm: Any = world, **kwargs: Any):
    """Factory matching ASE's ase.io.trajectory.Trajectory API.

    Returns TrajectoryReader in read mode, TrajectoryWriter otherwise.
    """
    if mode == 'r':
        return TrajectoryReader(filename)
    return TrajectoryWriter(filename, mode=mode, atoms=atoms, properties=properties, master=master, comm=comm, **kwargs)

