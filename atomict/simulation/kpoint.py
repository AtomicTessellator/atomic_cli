from atomict.api import get


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
