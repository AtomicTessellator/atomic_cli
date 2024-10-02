from atomict.api import delete, get, post


def create_simulation(
    project_id: str,
    control_file: str,
    geometry_file: str,
    action: str,
    name: str = None,
    description: str = None,
) -> dict:

    if action not in ["DRAFT", "LAUNCH"]:
        raise ValueError("Action must be 'DRAFT' or 'LAUNCH'")

    payload = {
        "project": project_id,
        "control_file": control_file,
        "geometry_file": geometry_file,
        "action": action,
        "name": name,
        "description": description,
    }

    result = post(
        "api/fhiaims-simulation/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def get_simulation(simulation_id: str):
    """
    Get a FHI aims simulation
    """
    result = get(f"api/fhiaims-simulation/{simulation_id}/")
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
        payload={"user_upload": user_upload_id, "simulation": fhi_simulation_id},
    )
    return result


def get_simulation_files(simulation_id: str):
    """
    Get the files associated with a FHI-aims simulation
    """
    result = get(f"api/fhiaims-simulation-file/?simulation_uuid={simulation_id}")
    return result
