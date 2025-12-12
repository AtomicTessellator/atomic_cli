from ase import Atoms
from ase.io import read

import logging
import os

# from atomict.io.msgpack import load_msgpack_trajectory
from atomict.io.formats.atraj import read_atraj
from atomict.io.formats.tess import read_tess
from atomict.io.fhiaims import read_aims_output
from atomict.io.utils import human_filesize
from atomict.simulation.mlrelax import get_mlrelax, get_mlrelax_files
from atomict.simulation.fhi_aims import get_simulation as fhi_get_simulation
from atomict.simulation.fhi_aims import get_simulation_files as fhi_get_simulation_files
from atomict.user.files import download_file
from atomict.user.workspace import download_workspace


def fetch_source_geometry(sim: dict, workbench_dir: str) -> Atoms:
    if sim.get("source_geometry"):

        extension = sim["source_geometry"]["orig_name"].split(".")[-1]

        download_file(sim["source_geometry"]["id"], workbench_dir + f"/geometry.{extension}")
        
        if extension == "atraj":
            atoms, _ = read_atraj(workbench_dir + f"/geometry.{extension}")
        elif extension == "tess":
            atoms, _ = read_tess(workbench_dir + f"/geometry.{extension}")
        else:
            atoms = read(workbench_dir + f"/geometry.{extension}")

        if isinstance(atoms, list):
            return atoms[-1]
        else:
            return atoms
    else:
        raise ValueError("No input geometry found")


def fetch_relaxed_geometry(sim: dict, workbench_dir: str) -> Atoms:

    """
    Fetch the relaxed geometry from the simulation
        sim can be any of these: FHIAimsSimulation, MLRelaxation, UserUpload
    
        returns: Atoms object
    """

    if sim.get("starting_structure"):
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

    elif sim.get("starting_structure_mlrelax"):
        
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
            atoms, _ = read_atraj(atraj_file)
        else:
            atoms = read(traj_file)
        
        if isinstance(atoms, list):
            return atoms[-1]
        else:
            return atoms

    elif sim.get("starting_structure_userupload"):
        logging.info(f"Previous UserUpload: {sim['starting_structure_userupload']['id']}")

        extension = sim["starting_structure_userupload"]["orig_name"].split(".")[-1]

        download_file(sim["starting_structure_userupload"]["id"], workbench_dir + f"/geometry.{extension}")
        
        if extension == "atraj":
            atoms, _ = read_atraj(workbench_dir + f"/geometry.{extension}")
        elif extension == "tess":
            atoms, _ = read_tess(workbench_dir + f"/geometry.{extension}")
        else:
            atoms = read(workbench_dir + f"/geometry.{extension}")

        if isinstance(atoms, list):
            return atoms[-1]
        else:
            return atoms
    else:
        raise ValueError("No input geometry found")
