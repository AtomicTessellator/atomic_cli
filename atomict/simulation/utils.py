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


logger = logging.getLogger(__name__)


def _read_geometry_file(filepath: str, extension: str) -> Atoms:
    if extension == "atraj":
        atoms, _ = read_atraj(filepath)
    elif extension == "tess":
        atoms, _ = read_tess(filepath)
    else:
        atoms = read(filepath)

    if isinstance(atoms, list):
        return atoms[-1]
    return atoms


def fetch_source_geometry(sim: dict, workbench_dir: str, *, api_root: str = None, token: str = None) -> Atoms:
    if sim.get("source_geometry"):
        extension = sim["source_geometry"]["orig_name"].split(".")[-1]
        filepath = workbench_dir + f"/geometry.{extension}"
        download_file(sim["source_geometry"]["id"], filepath, api_root=api_root, token=token)
        return _read_geometry_file(filepath, extension)
    else:
        raise ValueError("No associated input geometry found (simulation.source_geometry)")


def get_source_geometry(sim: dict, workbench_dir: str, *, api_root: str = None, token: str = None) -> Atoms:
    """Return the source geometry, reading from disk if already present."""
    if sim.get("source_geometry"):
        extension = sim["source_geometry"]["orig_name"].split(".")[-1]
        filepath = workbench_dir + f"/geometry.{extension}"
        if os.path.exists(filepath):
            logger.info(f"Source geometry already on disk, skipping download: {filepath}")
            return _read_geometry_file(filepath, extension)
        return fetch_source_geometry(sim, workbench_dir, api_root=api_root, token=token)
    else:
        raise ValueError("No associated input geometry found (simulation.source_geometry)")


def fetch_relaxed_geometry(sim: dict, workbench_dir: str, *, api_root: str = None, token: str = None) -> Atoms:

    """
    Fetch the relaxed geometry from the simulation
        sim can be any of these: FHIAimsSimulation, MLRelaxation, UserUpload
    
        returns: Atoms object
    """

    if sim.get("starting_structure"):
        previous_simulation = fhi_get_simulation(sim["starting_structure"]["id"], include_ht=True, api_root=api_root, token=token)
        logger.info(f"Previous simulation: {previous_simulation['id']}")
        files = fhi_get_simulation_files(previous_simulation["id"], api_root=api_root, token=token)

        total_size = 0
        for file in files["results"]:
            total_size += file["user_upload"]["size"]

        logger.info(
            f"Previous simulation: Downloading {len(files['results'])} files, Total size: {human_filesize(total_size)}"
        )

        prev_sim_dir = os.path.join(workbench_dir, "previous_simulation")
        os.makedirs(prev_sim_dir, exist_ok=True)
        download_workspace(files["results"], prev_sim_dir, api_root=api_root, token=token)
        atoms = read_aims_output(
            os.path.join(prev_sim_dir, f"{previous_simulation['id']}.out")
        )

        return atoms[-1]

    elif sim.get("starting_structure_mlrelax"):
        
        previous_mlrelax = get_mlrelax(sim["starting_structure_mlrelax"]["id"], include_ht=True, api_root=api_root, token=token)
        logger.info(f"Previous MLRelaxation: {previous_mlrelax['id']}")
        files = get_mlrelax_files(previous_mlrelax["id"], api_root=api_root, token=token)

        total_size = 0
        for file in files["results"]:
            total_size += file["user_upload"]["size"]

        logger.info(
            f"Previous MLRelaxation: Downloading {len(files['results'])} files, Total size: {human_filesize(total_size)}"
        )

        mlrelax_dir = os.path.join(workbench_dir, "previous_mlrelax")
        os.makedirs(mlrelax_dir, exist_ok=True)
        download_workspace(files["results"], mlrelax_dir, api_root=api_root, token=token)

        for ext in ("atraj", "tess", "traj"):
            candidate = os.path.join(mlrelax_dir, f"relax.{ext}")
            if os.path.exists(candidate):
                return _read_geometry_file(candidate, ext)

        raise FileNotFoundError(f"No relaxation output found in {mlrelax_dir}")

    elif sim.get("starting_structure_userupload"):
        logger.info(f"Previous UserUpload: {sim['starting_structure_userupload']['id']}")

        extension = sim["starting_structure_userupload"]["orig_name"].split(".")[-1]

        download_file(sim["starting_structure_userupload"]["id"], workbench_dir + f"/geometry.{extension}", api_root=api_root, token=token)
        
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
