try:
    from ase import Atoms
    from ase.io import read
except ImportError:
    raise ImportError(
        "The 'ase' package is required for JSON operations with Atoms objects. "
        "To install the optional dependencies such as ase, spglib, pymatgen, use `pip install atomict[utils]`"
    )

import logging
import os

from atomict.io.msgpack import load_msgpack_trajectory
from atomict.io.fhiaims import read_aims_output
from atomict.io.utils import human_filesize
from atomict.simulation.mlrelax import get_mlrelax, get_mlrelax_files
from atomict.simulation.fhi_aims import get_simulation as fhi_get_simulation
from atomict.simulation.fhi_aims import get_simulation_files as fhi_get_simulation_files
from atomict.user.files import download_file
from atomict.user.workspace import download_workspace


def fetch_relaxed_geometry(sim: dict, workbench_dir: str) -> Atoms:

    """
    Fetch the relaxed geometry from the simulation
        sim can be any of these: FHIAimsSimulation, MLRelaxation, UserUpload
    
        returns: Atoms object
    """

    if sim["starting_structure"]:
        previous_simulation = fhi_get_simulation(sim["starting_structure"]["id"], include_ht=True)
        logging.info(f"Previous simulation: {previous_simulation['id']}")
        files = fhi_get_simulation_files(previous_simulation["id"])

        total_size = 0
        for file in files["results"]:
            total_size += file["user_upload"]["size"]

        logging.info(
            f"Previous simulation: Downloading {len(files['results'])} files, Total size: {human_filesize(total_size)}"
        )

        prev_sim_dir = os.path.join(workbench_dir, "previous_simulation")
        os.makedirs(prev_sim_dir, exist_ok=True)
        download_workspace(files["results"], prev_sim_dir)
        atoms = read_aims_output(
            os.path.join(prev_sim_dir, f"{previous_simulation['id']}.out")
        )

        return atoms[-1]

    elif sim["starting_structure_mlrelax"]:
        
        previous_mlrelax = get_mlrelax(sim["starting_structure_mlrelax"]["id"], include_ht=True)
        logging.info(f"Previous MLRelaxation: {previous_mlrelax['id']}")
        files = get_mlrelax_files(previous_mlrelax["id"])

        total_size = 0
        for file in files["results"]:
            total_size += file["user_upload"]["size"]

        logging.info(
            f"Previous MLRelaxation: Downloading {len(files['results'])} files, Total size: {human_filesize(total_size)}"
        )

        mlrelax_dir = os.path.join(workbench_dir, "previous_mlrelax")
        os.makedirs(mlrelax_dir, exist_ok=True)
        download_workspace(files["results"], mlrelax_dir)

        traj_file = os.path.join(mlrelax_dir, "relax.traj")
        atraj_file = os.path.join(mlrelax_dir, "relax.atraj")

        if os.path.exists(atraj_file):
            atoms, _ = load_msgpack_trajectory(atraj_file)
        else:
            atoms = read(traj_file)
        
        if isinstance(atoms, list):
            return atoms[-1]
        else:
            return atoms

    elif sim["starting_structure_userupload"]:
        logging.info(f"Previous UserUpload: {sim['starting_structure_userupload']['uuid']}")

        extension = sim["starting_structure_userupload"]["orig_name"].split(".")[-1]

        download_file(sim["starting_structure_userupload"]["uuid"], workbench_dir + f"/relaxed.{extension}")
        atoms = read(workbench_dir + f"/relaxed.{extension}")

        if isinstance(atoms, list):
            return atoms[-1]
        else:
            return atoms
    else:
        raise ValueError("No relaxed structure simulation found")
