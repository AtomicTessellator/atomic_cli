import os
import logging

try:
    from ase import Atoms
    from ase.io import read
    from ase.io.formats import UnknownFileTypeError
except ImportError:
    raise ImportError(
        "The 'ase' package is required for JSON operations with Atoms objects. "
        "To install the optional dependencies such as ase, spglib, pymatgen, use `pip install atomict[utils]`"
    )

from atomict.io.fhiaims import read_aims_output
from atomict.utils.path_security import ensure_safe_read_path, PathTraversalError


def read_final_geometry(workspace_dir: str, simulation_id: str):
    """Read the starting structure from the workspace directory

    Args:
        workspace_dir (str): The workspace directory
    """
    next_step_path = f"{workspace_dir}/starting_structure/geometry.in.next_step"
    if os.path.exists(next_step_path):
        try:
            return read(next_step_path, foramt='aims')
        except UnknownFileTypeError:
            logging.warning(f"Could not read {next_step_path}, falling back to output file")

    return read_aims_output(
        f"{workspace_dir}/starting_structure/{simulation_id}.out"
    )[-1]


def fhi_to_ase(geometry_file):
    """
    Convert FHI-aims geometry.in to ASE atoms object.
    Handles both periodic and non-periodic systems.
    
    Args:
        geometry_file (str): Path to FHI-aims geometry.in file (validated for security)
        
    Returns:
        ase.Atoms: ASE atoms object with appropriate periodicity
        
    Raises:
        PathTraversalError: If the path contains traversal attempts
        FileNotFoundError: If the file doesn't exist
    """
    # Validate the file path for safe reading
    safe_path = ensure_safe_read_path(geometry_file)
    
    lattice_vectors = []
    positions = []
    symbols = []
    
    logging.debug(f"Reading FHI-aims geometry from: {safe_path}")
    with open(safe_path, 'r') as f:
        for line in f:
            if line.startswith('lattice_vector'):
                vec = [float(x) for x in line.split()[1:4]]
                lattice_vectors.append(vec)
            elif line.startswith('atom'):
                parts = line.split()
                pos = [float(x) for x in parts[1:4]]
                symbol = parts[4]
                positions.append(pos)
                symbols.append(symbol)

    # Check if the system is periodic (has lattice vectors)
    if lattice_vectors:
        atoms = Atoms(symbols=symbols, 
                    positions=positions, 
                    cell=lattice_vectors, 
                    pbc=True)
    else:
        atoms = Atoms(symbols=symbols, 
                    positions=positions, 
                    pbc=False)

    return atoms
