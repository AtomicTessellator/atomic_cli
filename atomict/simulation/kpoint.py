from atomict.api import get, post


def get_kpoint_exploration(simulation_id: str):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-exploration/{simulation_id}/")
    return result


def get_kpoint_simulation(simulation_id: str):
    """
    Get KPoint simulation
    """
    result = get(f"api/kpoint-simulation/{simulation_id}/")
    return result


def create_kpoint_simulation(exploration_id: str, simulation_id: str, k_points: list[float]):
    """
    Create KPoint simulation
    """
    result = post(
        f"api/kpoint-simulation/",
        data={"exploration_id": exploration_id, "simulation_id": simulation_id, "k_points": k_points},
    )
    return result
