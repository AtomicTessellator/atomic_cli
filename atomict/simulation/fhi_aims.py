from atomict.api import delete, get, post
from atomict.infra.distwork.task import SimulationAction


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


def get_simulation(simulation_id: str, **params):
    """
    Get a FHI aims simulation
    """
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/fhiaims-simulation/{simulation_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    result = get(url)
    return result


def delete_simulation(simulation_id):
    """
    Delete a FHI aims simulation
    """
    result = delete(f"api/fhiaims-simulation/{simulation_id}/")
    return result


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
    result = get(f"api/fhiaims-simulation-file/?simulation__id={simulation_id}")
    return result
