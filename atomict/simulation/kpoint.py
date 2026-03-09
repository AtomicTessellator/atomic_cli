from atomict.api import get, post, patch
from atomict.resource_helpers import delete_resource, list_resources, retrieve_resource


def get_kpoint_exploration(simulation_id: str, **params):
    """Get K-point exploration details.
    
    Args:
        simulation_id (str): The simulation identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/kpoint-exploration", simulation_id, **params)


def list_kpoint_explorations(**params):
    """List K-point explorations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/kpoint-exploration", **params)


def create_kpoint_exploration(payload: dict[str, object]):
    """Create a new K-point exploration.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/kpoint-exploration/", payload=payload)


def update_kpoint_exploration(exploration_id: str, fields: dict):
    """Update K-point exploration.
    
    Args:
        exploration_id (str): The exploration identifier.
        fields (dict): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    result = patch(f"api/kpoint-exploration/{exploration_id}/", payload=fields)
    return result


def delete_kpoint_exploration(exploration_id: str):
    """Delete K-point exploration.
    
    Args:
        exploration_id (str): The exploration identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/kpoint-exploration", exploration_id)


def get_kpoint_simulation_list(exploration_id: str):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-simulation/?exploration__id={exploration_id}")
    return result


def list_kpoint_simulations(**params):
    """List K-point simulations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/kpoint-simulation", **params)


def get_kpoint_simulation(simulation_id: str, **params):
    """Get K-point simulation details.
    
    Args:
        simulation_id (str): The simulation identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/kpoint-simulation", simulation_id, **params)


def create_kpoint_simulation(
    exploration_id: str, simulation_id: str, k_points: list[float]
):
    """Create a new K-point simulation.
    
    Args:
        exploration_id (str): The exploration identifier.
        simulation_id (str): The simulation identifier.
        k_points (list[float]): The K-point values for the simulation.
    
    Returns:
        dict: The API response payload.
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
    """Update K-point simulation.
    
    Args:
        simulation_id (str): The simulation identifier.
        fields (dict): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    result = patch(f"api/kpoint-simulation/{simulation_id}/", payload=fields)
    return result


def delete_kpoint_simulation(simulation_id: str):
    """Delete K-point simulation.
    
    Args:
        simulation_id (str): The simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/kpoint-simulation", simulation_id)


def get_kpoint_analysis(analysis_id: str, **params):
    """Get K-point analysis details.
    
    Args:
        analysis_id (str): The analysis identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/kpoint-analysis", analysis_id, **params)


def list_kpoint_analyses(**params):
    """List K-point analyses.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/kpoint-analysis", **params)


def create_kpoint_analysis(payload: dict[str, object]):
    """Create a new K-point analysis.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/kpoint-analysis/", payload=payload)


def update_kpoint_analysis(analysis_id: str, fields: dict):
    """Update K-point analysis.
    
    Args:
        analysis_id (str): The analysis identifier.
        fields (dict): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    result = patch(f"api/kpoint-analysis/{analysis_id}/", payload=fields)
    return result


def delete_kpoint_analysis(analysis_id: str):
    """Delete K-point analysis.
    
    Args:
        analysis_id (str): The analysis identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/kpoint-analysis", analysis_id)
