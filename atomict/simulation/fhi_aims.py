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
