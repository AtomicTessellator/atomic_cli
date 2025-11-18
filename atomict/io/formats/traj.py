from typing import Union, List, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    import ase


def write_traj(atoms: Union['ase.Atoms', List['ase.Atoms']], filename: str, metadata: Dict = None, properties: List[str] = None):
    """Save atoms to an ASE trajectory file.
    
    Parameters:
    atoms: Single Atoms object or list of Atoms objects to write
    filename: Output trajectory filename
    metadata: Optional metadata to add to frames
    properties: List of calculator properties to save (e.g., ['energy', 'forces', 'stress'])
                If None, all available properties are saved
    """
    
    from ase.io.trajectory import TrajectoryWriter
    from ase import Atoms
    
    # Single atoms case - convert to list
    if isinstance(atoms, Atoms):
        atoms_list = [atoms]
    else:
        atoms_list = atoms
    
    # Use TrajectoryWriter for better control over properties
    with TrajectoryWriter(filename, mode='w', properties=properties) as traj:
        for atoms_frame in atoms_list:
            if metadata:
                # Only copy if we need to add metadata to avoid modifying original
                frame_copy = atoms_frame.copy()
                frame_copy.info.update(metadata)
                traj.write(frame_copy)
            else:
                # Write directly without copying if no metadata
                traj.write(atoms_frame)


def read_traj(filename: str) -> tuple[List['ase.Atoms'], Dict]:
    """Load atoms from an ASE trajectory file."""
    
    from ase.io.trajectory import TrajectoryReader
    
    atoms_list = []
    metadata = {}
    
    # Read trajectory using TrajectoryReader for better control
    with TrajectoryReader(filename) as traj:
        for atoms_frame in traj:
            atoms_list.append(atoms_frame)
            # Extract metadata from the first frame's info
            if not metadata and atoms_frame.info:
                metadata = dict(atoms_frame.info)
    
    return atoms_list, metadata
