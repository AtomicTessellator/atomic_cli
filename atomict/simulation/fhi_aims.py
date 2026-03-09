from atomict.api import post
from atomict.infra.distwork.task import SimulationAction
from atomict.resource_helpers import (
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def create_simulation(
    project_id: str,
    control_file: str,
    geometry_file: str,
    action: SimulationAction,
    name: str = None,
    description: str = None,
    extra_simulation_kwargs: dict = None,
) -> dict:

    """Create a new FHI-aims simulation.
    
    Args:
        project_id (str): The project identifier.
        control_file (str): The FHI-aims control file contents.
        geometry_file (str): The FHI-aims geometry file contents.
        action (SimulationAction): The action to perform for the simulation.
        name (str | None): The resource name.
        description (str | None): The resource description.
        extra_simulation_kwargs (dict | None): Additional simulation keyword arguments to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    if action not in [SimulationAction.SAVE_DRAFT, SimulationAction.LAUNCH]:
        raise ValueError("Action must be 'SimulationAction.SAVE_DRAFT' or 'SimulationAction.LAUNCH'")

    payload = {
        "project_id": project_id,
        "control_file": control_file,
        "geometry_file": geometry_file,
        "action": action.value,
        "name": name,
        "description": description,
    }

    if extra_simulation_kwargs:
        payload.update(extra_simulation_kwargs)

    result = post(
        "api/fhiaims-simulation/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def create_fhiaims_simulation(*args, **kwargs):
    return create_simulation(*args, **kwargs)


def get_simulation(simulation_id: str, **params):
    """Get FHI-aims simulation details.
    
    Args:
        simulation_id (str): The simulation identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/fhiaims-simulation", simulation_id, **params)


def get_fhiaims_simulation(simulation_id: str, **params):
    return get_simulation(simulation_id, **params)


def list_simulations(**params):
    """List FHI-aims simulations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/fhiaims-simulation", **params)


def list_fhiaims_simulations(**params):
    return list_simulations(**params)


def update_simulation(simulation_id: str, fields: dict[str, object]):
    """Update FHI-aims simulation.
    
    Args:
        simulation_id (str): The simulation identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/fhiaims-simulation", simulation_id, fields)


def update_fhiaims_simulation(simulation_id: str, fields: dict[str, object]):
    return update_simulation(simulation_id, fields)


def delete_simulation(simulation_id):
    """Delete FHI-aims simulation.
    
    Args:
        simulation_id (str): The simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/fhiaims-simulation", simulation_id)


def delete_fhiaims_simulation(simulation_id: str):
    return delete_simulation(simulation_id)


def associate_user_upload_with_fhiaims_simulation(
    user_upload_id: str, fhi_simulation_id: str
):
    """
    Associate a user upload with a FHI-aims simulation
    """
    result = post(
        "api/fhiaims-simulation-file/",
        payload={"user_upload_id": user_upload_id, "simulation_id": fhi_simulation_id},
    )
    return result


def get_simulation_files(simulation_id: str):
    """
    Get the files associated with a FHI-aims simulation
    """
    return list_resources("api/fhiaims-simulation-file", simulation__id=simulation_id)


def get_simulation_file(file_id: str, **params):
    """Get FHI-aims simulation file details.
    
    Args:
        file_id (str): The file identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/fhiaims-simulation-file", file_id, **params)


def get_fhiaims_simulation_file(file_id: str, **params):
    return get_simulation_file(file_id, **params)


def list_simulation_files(**params):
    """List FHI-aims simulation files.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/fhiaims-simulation-file", **params)


def list_fhiaims_simulation_files(**params):
    return list_simulation_files(**params)


def create_simulation_file(payload: dict[str, object]):
    """Create a new FHI-aims simulation file.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/fhiaims-simulation-file/", payload=payload)


def update_simulation_file(file_id: str, fields: dict[str, object]):
    """Update FHI-aims simulation file.
    
    Args:
        file_id (str): The file identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/fhiaims-simulation-file", file_id, fields)


def delete_simulation_file(file_id: str):
    """Delete FHI-aims simulation file.
    
    Args:
        file_id (str): The file identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/fhiaims-simulation-file", file_id)
