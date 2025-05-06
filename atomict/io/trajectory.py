"""Trajectory implementation using msgpack storage"""
from typing import Union, List, Dict, Any, Optional, Tuple
import os
import contextlib
import io
import warnings
import numpy as np

from ase import __version__
from ase.atoms import Atoms
from ase.calculators.calculator import PropertyNotImplementedError
from ase.calculators.singlepoint import SinglePointCalculator, all_properties
from ase.io.formats import is_compressed
from ase.constraints import dict2constraint
from ase.io.jsonio import decode, encode
from ase.parallel import world

__all__ = ['Trajectory']


def Trajectory(filename, mode='r', atoms=None, properties=None, master=None, comm=world, flush_interval=1):
    """A Trajectory can be created in read, write or append mode.

    Parameters:

    filename: str
        The name of the file.  Traditionally ends in .traj, but .mpk is recommended for msgpack.
    mode: str
        The mode.  'r' is read mode, the file should already exist, and
        no atoms argument should be specified.
        'w' is write mode.  The atoms argument specifies the Atoms
        object to be written to the file, if not given it must instead
        be given as an argument to the write() method.
        'a' is append mode.  It acts as write mode, except that
        data is appended to a preexisting file.
    atoms: Atoms object
        The Atoms object to be written in write or append mode.
    properties: list of str
        If specified, these calculator properties are saved in the
        trajectory.  If not specified, all supported quantities are
        saved.  Possible values: energy, forces, stress, dipole,
        charges, magmom and magmoms.
    master: bool
        Controls which process does the actual writing. The
        default is that process number 0 does this.  If this
        argument is given, processes where it is True will write.
    comm: Communicator object
        Communicator to handle parallel file reading and writing.
    flush_interval: int
        Controls how often the trajectory is written to disk. By default,
        the trajectory is written to disk after every frame (flush_interval=1).
        For production runs with many frames, a higher value like 10 or 100
        can significantly improve performance by reducing I/O operations.

    The atoms, properties and master arguments are ignored in read mode.
    """
    if mode == 'r':
        return TrajectoryReader(filename)
    return TrajectoryWriter(filename, mode, atoms, properties, master=master, comm=comm, flush_interval=flush_interval)


class TrajectoryWriter:
    """Writes Atoms objects to a .mpk file using msgpack."""

    def __init__(self, filename, mode='w', atoms=None, properties=None, master=None, comm=world, flush_interval=1):
        """A Trajectory writer, in write or append mode.

        Parameters:

        filename: str
            The name of the file. .mpk extension is recommended for msgpack.
        mode: str
            The mode.  'w' is write mode.  The atoms argument specifies the Atoms
            object to be written to the file, if not given it must instead
            be given as an argument to the write() method.
            'a' is append mode.  It acts as write mode, except that
            data is appended to a preexisting file.
        atoms: Atoms object
            The Atoms object to be written in write or append mode.
        properties: list of str
            If specified, these calculator properties are saved in the
            trajectory.  If not specified, all supported quantities are
            saved.  Possible values: energy, forces, stress, dipole,
            charges, magmom and magmoms.
        master: bool
            Controls which process does the actual writing. The
            default is that process number 0 does this.  If this
            argument is given, processes where it is True will write.
        comm: MPI communicator
            MPI communicator for this trajectory writer, by default world.
        flush_interval: int
            Controls how often the trajectory is written to disk. By default,
            the trajectory is written to disk after every frame (flush_interval=1).
            For production runs with many frames, a higher value like 10 or 100
            can significantly improve performance by reducing I/O operations.
        """
        try:
            import msgpack
            import msgpack_numpy as m
        except ImportError:
            raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")
        
        # Enable numpy array serialization
        m.patch()
        
        if master is None:
            master = comm.rank == 0

        self.filename = filename
        self.mode = mode
        self.atoms = atoms
        self.properties = properties
        self.master = master
        self.comm = comm
        self.description = {}
        self.flush_interval = max(1, flush_interval)  # Ensure at least 1
        
        # Track if this is first write to ensure it's always flushed
        self._is_first_write = True
        self._frames = []
        self._frames_since_last_flush = 0
        self._total_frames_added = 0  # Track total frames added for flush interval
        
        # Store metadata in separate dict to ensure it's preserved
        self._metadata = {
            "description": {},
            "ase_version": __version__
        }
        
        self._open(filename, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.close()

    def set_description(self, description):
        """Set the description metadata for the trajectory."""
        self.description.update(description)
        self._metadata["description"].update(description)

    def _open(self, filename, mode):
        if mode not in 'aw':
            raise ValueError('mode must be "w" or "a".')
            
        if mode == 'a' and os.path.exists(self.filename) and self.master:
            # Load existing data in append mode
            with TrajectoryReader(self.filename) as reader:
                self._frames = [atoms for atoms in reader]
                # Copy description if available
                if reader.description:
                    self.description.update(reader.description)
                    self._metadata["description"].update(reader.description)
                self._is_first_write = False
                self._total_frames_added = len(self._frames)
        else:
            self._frames = []
            self._is_first_write = True
            self._total_frames_added = 0

    def write(self, atoms=None, **kwargs):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.

        Use keyword arguments to add extra properties::

            writer.write(atoms, energy=117, dipole=[0, 0, 1.0])
        """
        if atoms is None:
            atoms = self.atoms

        for image in atoms.iterimages():
            self._write_atoms(image, **kwargs)

    def _write_atoms(self, atoms, **kwargs):
        if not self.master:
            return
            
        # Make a copy of the atoms to save
        atoms_copy = atoms.copy()
        
        # Handle calculator data
        calc = atoms.calc
        
        if calc is None and len(kwargs) > 0:
            calc = SinglePointCalculator(atoms)
            
        if calc is not None:
            # Store calculator results in atoms info
            calc_data = {'name': calc.name}
            
            if hasattr(calc, 'todict'):
                calc_data['parameters'] = calc.todict()
                
            for prop in all_properties:
                if prop in kwargs:
                    x = kwargs[prop]
                elif self.properties is not None:
                    if prop in self.properties:
                        try:
                            x = calc.get_property(prop, atoms)
                        except (PropertyNotImplementedError, KeyError):
                            x = None
                    else:
                        x = None
                else:
                    try:
                        x = calc.get_property(prop, atoms, allow_calculation=False)
                    except (PropertyNotImplementedError, KeyError):
                        x = None
                        
                if x is not None:
                    if prop in ['stress', 'dipole']:
                        x = x.tolist()
                    calc_data[prop] = x
                    
            atoms_copy.info['_calc_data'] = calc_data
            
        # Store constraints as encoded data
        if atoms_copy.constraints:
            atoms_copy.info['_constraints'] = encode(atoms_copy.constraints)
            
        # Always store description and ASE version in atoms info as well
        atoms_copy.info['_traj_description'] = self.description
        atoms_copy.info['_ase_version'] = __version__
            
        self._frames.append(atoms_copy)
        self._frames_since_last_flush += 1
        self._total_frames_added += 1
        
        # Save to disk if:
        # 1. First write (always flush first frame)
        # 2. We've accumulated enough frames based on flush interval
        if self._is_first_write or self._total_frames_added % self.flush_interval == 0:
            self._save()
            self._frames_since_last_flush = 0
            self._is_first_write = False

    def _save(self):
        """Save all frames to the file."""
        if not self.master or not self._frames:
            return
            
        from .msgpack import save_msgpack_trajectory
        
        # Save with metadata
        save_msgpack_trajectory(self._frames, self.filename, metadata=self._metadata)

    def flush(self):
        """Flush the trajectory to disk."""
        if self.master and self._frames_since_last_flush > 0:
            self._save()
            self._frames_since_last_flush = 0

    def close(self):
        """Close the trajectory file."""
        if self.master and self._frames_since_last_flush > 0:
            self._save()
        self._frames_since_last_flush = 0

    def __len__(self):
        return self.comm.sum_scalar(len(self._frames))


class TrajectoryReader:
    """Reads Atoms objects from a msgpack trajectory file."""

    def __init__(self, filename):
        """A Trajectory in read mode.

        The filename traditionally ends in .traj, but .mpk is recommended for msgpack.
        """
        try:
            import msgpack
            import msgpack_numpy as m
        except ImportError:
            raise ImportError("You need to install with `pip install atomict[tools]` to use msgpack I/O")
            
        # Enable numpy array deserialization
        m.patch()
        
        self.filename = filename
        self._frames = None
        self.description = None
        self.ase_version = None
        
        self._open(filename)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.close()

    def _open(self, filename):
        from .msgpack import load_msgpack_trajectory
        
        # Load trajectory with metadata
        result = load_msgpack_trajectory(filename)
        
        # Handle various ways the msgpack loader might return data
        if isinstance(result, tuple) and len(result) == 2:
            atoms_list, metadata = result
        elif isinstance(result, dict) and 'frames' in result and 'metadata' in result:
            atoms_list = result['frames']
            metadata = result['metadata']
        else:
            # Fall back for backward compatibility
            atoms_list = result
            metadata = {}
        
        # Make sure we have a list
        if not isinstance(atoms_list, list):
            atoms_list = [atoms_list]
            
        self._frames = atoms_list
        
        # Extract description from metadata first, then fallback to atoms info
        if metadata and 'description' in metadata and metadata['description']:
            self.description = metadata['description']
        elif len(self._frames) > 0:
            # Try various description keys in atoms info
            for key in ['_traj_description', '_description']:
                if key in self._frames[0].info:
                    self.description = self._frames[0].info[key]
                    break
            
        # Extract version
        if metadata and 'ase_version' in metadata:
            self.ase_version = metadata['ase_version']
        elif len(self._frames) > 0 and '_ase_version' in self._frames[0].info:
            self.ase_version = self._frames[0].info['_ase_version']
        else:
            self.ase_version = 'unknown'

    def close(self):
        """Close the trajectory file."""
        # Nothing to do for msgpack files
        pass

    def __getitem__(self, i=-1):
        if isinstance(i, slice):
            return SlicedTrajectory(self, i)
            
        atoms = self._frames[i].copy()
        
        # Process calculator data if present
        if '_calc_data' in atoms.info:
            calc_data = atoms.info.pop('_calc_data')
            
            results = {}
            implemented_properties = []
            
            for prop in all_properties:
                if prop in calc_data:
                    results[prop] = calc_data[prop]
                    implemented_properties.append(prop)
                    
            calc = SinglePointCalculator(atoms, **results)
            calc.name = calc_data.get('name', 'unknown')
            calc.implemented_properties = implemented_properties
            
            if 'parameters' in calc_data:
                calc.parameters.update(calc_data['parameters'])
                
            atoms.calc = calc
            
        # Process constraints if present
        if '_constraints' in atoms.info:
            constraints_data = atoms.info.pop('_constraints')
            atoms.constraints = [dict2constraint(d) for d in decode(constraints_data)]
            
        # Clean up any remaining metadata keys
        for key in ['_traj_description', '_description', '_ase_version']:
            if key in atoms.info:
                atoms.info.pop(key)
            
        return atoms

    def __len__(self):
        return len(self._frames)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


class SlicedTrajectory:
    """Wrapper to return a slice from a trajectory without loading
    from disk. Initialize with a trajectory (in read mode) and the
    desired slice object."""

    def __init__(self, trajectory, sliced):
        self.trajectory = trajectory
        self.map = range(len(self.trajectory))[sliced]

    def __getitem__(self, i):
        if isinstance(i, slice):
            # Map directly to the original traj, not recursively.
            traj = SlicedTrajectory(self.trajectory, slice(0, None))
            traj.map = self.map[i]
            return traj
        return self.trajectory[self.map[i]]

    def __len__(self):
        return len(self.map)


def read_traj(fd, index):
    """Read msgpack trajectory for ase.io.read()."""
    trj = TrajectoryReader(fd)
    for i in range(*index.indices(len(trj))):
        yield trj[i]


@contextlib.contextmanager
def defer_compression(fd):
    """Defer the file compression until all the configurations are read."""
    # We do this because the trajectory and compressed-file
    # internals do not play well together.
    # Be advised not to defer compression of very long trajectories
    # as they use a lot of memory.
    if is_compressed(fd):
        with io.BytesIO() as bytes_io:
            try:
                # write the uncompressed data to the buffer
                yield bytes_io
            finally:
                # write the buffered data to the compressed file
                bytes_io.seek(0)
                fd.write(bytes_io.read())
    else:
        yield fd


def write_traj(fd, images):
    """Write image(s) to msgpack trajectory."""
    if isinstance(images, Atoms):
        images = [images]
    with defer_compression(fd) as fd_uncompressed:
        trj = TrajectoryWriter(fd_uncompressed)
        for atoms in images:
            trj.write(atoms)
