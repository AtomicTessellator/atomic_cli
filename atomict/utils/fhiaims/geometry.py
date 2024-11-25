from atomict.io.fhiaims import read_aims_output
from ase.io import read
from ase.io.formats import UnknownFileTypeError
import os
import logging


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
