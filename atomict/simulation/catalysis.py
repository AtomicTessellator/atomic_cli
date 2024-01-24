from atomict.api import get


def get_exploration(exploration_id: str):
    return get(f"api/catalysis-exploration/{exploration_id}/")


def get_simulation(simulation_id: str):
    return get(f"api/catalysis-simulation/{simulation_id}/")
