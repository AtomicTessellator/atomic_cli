from atomict.api import delete, get, post


def create_simulation(
    project_id: str,
    input_file: str,
    action: str,
    name: str = None,
    description: str = None,
    **kwargs,
) -> dict:

    if action not in ["DRAFT", "LAUNCH"]:
        raise ValueError("Action must be 'DRAFT' or 'LAUNCH'")

    payload = {
        "project": project_id,
        "input_file": input_file,
        "action": action,
        "name": name,
        "description": description,
    }

    payload.update(kwargs)

    result = post(
        "api/qe-simulation/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def get_simulation_by_name(name: str) -> dict:
    result = get(f"api/qe-simulation/?name={name}")
    return result["results"][0]


def get_simulation(simulation_id: str):
    """
    Get a Quantum Espresso simulation
    """
    result = get(f"api/qe-simulation/{simulation_id}/")
    return result


def delete_simulation(simulation_id):
    """
    Delete a Quantum Espresso simulation
    """
    result = delete(f"api/qe-simulation/{simulation_id}/")
    return result


def associate_user_upload_with_qe_simulation(
    user_upload_id: str, qe_simulation_id: str
):
    """
    Associate a user upload with a Quantum Espresso simulation
    """
    result = post(
        "api/qe-simulation-file/",
        payload={"user_upload_id": user_upload_id, "simulation_id": qe_simulation_id},
    )
    return result


def get_simulation_files(simulation_id: str):
    """
    Get the files associated with a Quantum Espresso simulation
    """
    result = get(f"api/qe-simulation-file/?simulation_id={simulation_id}")
    return result
