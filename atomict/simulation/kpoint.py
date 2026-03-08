from atomict.api import get, post, patch
from atomict.resource_helpers import delete_resource, list_resources, retrieve_resource


def get_kpoint_exploration(simulation_id: str, **params):
    """
    Get kpoints for a simulation
    """
    return retrieve_resource("api/kpoint-exploration", simulation_id, **params)


def list_kpoint_explorations(**params):
    return list_resources("api/kpoint-exploration", **params)


def create_kpoint_exploration(payload: dict[str, object]):
    return post("api/kpoint-exploration/", payload=payload)


def update_kpoint_exploration(exploration_id: str, fields: dict):
    """
    Update KPoint exploration
    """
    result = patch(f"api/kpoint-exploration/{exploration_id}/", payload=fields)
    return result


def delete_kpoint_exploration(exploration_id: str):
    return delete_resource("api/kpoint-exploration", exploration_id)


def get_kpoint_simulation_list(exploration_id: str):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-simulation/?exploration__id={exploration_id}")
    return result


def list_kpoint_simulations(**params):
    return list_resources("api/kpoint-simulation", **params)


def get_kpoint_simulation(simulation_id: str, **params):
    """
    Get KPoint simulation
    """
    return retrieve_resource("api/kpoint-simulation", simulation_id, **params)


def create_kpoint_simulation(
    exploration_id: str, simulation_id: str, k_points: list[float]
):
    """
    Create KPoint simulation
    """
    result = post(
        "api/kpoint-simulation/",
        payload={
            "exploration_id": exploration_id,
            "simulation_id": simulation_id,
            "k_points": k_points,
        },
    )
    return result


def update_kpoint_simulation(simulation_id: str, fields: dict):
    """
    Update KPoint simulation
    """
    result = patch(f"api/kpoint-simulation/{simulation_id}/", payload=fields)
    return result


def delete_kpoint_simulation(simulation_id: str):
    return delete_resource("api/kpoint-simulation", simulation_id)


def get_kpoint_analysis(analysis_id: str, **params):
    """
    Get KPoint analysis
    """
    return retrieve_resource("api/kpoint-analysis", analysis_id, **params)


def list_kpoint_analyses(**params):
    return list_resources("api/kpoint-analysis", **params)


def create_kpoint_analysis(payload: dict[str, object]):
    return post("api/kpoint-analysis/", payload=payload)


def update_kpoint_analysis(analysis_id: str, fields: dict):
    """
    Update KPoint analysis
    """
    result = patch(f"api/kpoint-analysis/{analysis_id}/", payload=fields)
    return result


def delete_kpoint_analysis(analysis_id: str):
    return delete_resource("api/kpoint-analysis", analysis_id)
