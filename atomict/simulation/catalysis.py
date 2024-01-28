from atomict.api import get, post


def get_exploration(exploration_id: str):
    return get(f"api/catalysis-exploration/{exploration_id}/")


def get_simulation(simulation_id: str):
    return get(f"api/catalysis-simulation/{simulation_id}/")


def upload_simulation_results(results, traj_file_path, log_file_path):
    url = "simulation/catalysis/mlacc/results_upload/"

    files = {
        "traj": ("simulation.traj", open(traj_file_path, "rb")),
        "logs": ("simulation.log", open(log_file_path, "rb")),
    }

    return post(url, results, files=files)
