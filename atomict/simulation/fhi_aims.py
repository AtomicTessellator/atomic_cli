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
    """
    Get a FHI aims simulation
    """
    return retrieve_resource("api/fhiaims-simulation", simulation_id, **params)


def get_fhiaims_simulation(simulation_id: str, **params):
    return get_simulation(simulation_id, **params)


def list_simulations(**params):
    return list_resources("api/fhiaims-simulation", **params)


def list_fhiaims_simulations(**params):
    return list_simulations(**params)


def update_simulation(simulation_id: str, fields: dict[str, object]):
    return update_resource("api/fhiaims-simulation", simulation_id, fields)


def update_fhiaims_simulation(simulation_id: str, fields: dict[str, object]):
    return update_simulation(simulation_id, fields)


def delete_simulation(simulation_id):
    """
    Delete a FHI aims simulation
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
    return retrieve_resource("api/fhiaims-simulation-file", file_id, **params)


def get_fhiaims_simulation_file(file_id: str, **params):
    return get_simulation_file(file_id, **params)


def list_simulation_files(**params):
    return list_resources("api/fhiaims-simulation-file", **params)


def list_fhiaims_simulation_files(**params):
    return list_simulation_files(**params)


def create_simulation_file(payload: dict[str, object]):
    return post("api/fhiaims-simulation-file/", payload=payload)


def update_simulation_file(file_id: str, fields: dict[str, object]):
    return update_resource("api/fhiaims-simulation-file", file_id, fields)


def delete_simulation_file(file_id: str):
    return delete_resource("api/fhiaims-simulation-file", file_id)
