from atomict.io.fhiaims import read_aims_output
from ase.io import read
import os


def read_final_geometry(workspace_dir: str, simulation_id: str):
    """Read the starting structure from the workspace directory

    Args:
        workspace_dir (str): The workspace directory
    """
    if os.path.exists(f"{workspace_dir}/starting_structure/geometry.in.next_step"):
        return read(f"{workspace_dir}/starting_structure/geometry.in.next_step")
    else:
        return read_aims_output(
            f"{workspace_dir}/starting_structure/{simulation_id}.out"
        )[-1]
