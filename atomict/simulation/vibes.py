from typing import Any, Dict, Optional

from atomict.api import delete, get, post


def create_vibes_simulation(
    project_id: str,
    starting_structure_id: str,
    action: str = "DRAFT",
    name: Optional[str] = None,
    description: Optional[str] = None,
    calculator_parameters: Optional[Dict[str, Any]] = None,
    calculator_kpoints: Optional[Dict[str, Any]] = None,
    calculator_basis_set: Optional[Dict[str, Any]] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Any:
    """
    Create a VIBES phonon simulation

    Args:
        project_id: str - The ID of the project
        starting_structure_id: str - The ID of the FHI-aims starting structure
        action: str - The action to perform ('DRAFT' or 'LAUNCH')
        name: str - Optional name for the simulation
        description: str - Optional description for the simulation
        calculator_parameters: dict - Calculator parameters (e.g., {'xc': 'pw-lda'})
        calculator_kpoints: dict - K-point settings (e.g., {'density': 3.5})
        calculator_basis_set: dict - Basis set settings (e.g., {'default': 'light'})
        extra_kwargs: dict - Additional parameters for launch (e.g., cluster config)

    Returns:
        dict: The created VIBES simulation data

    Raises:
        ValueError: If action is not 'DRAFT' or 'LAUNCH'
    """
    if action not in ["DRAFT", "LAUNCH"]:
        raise ValueError("Action must be 'DRAFT' or 'LAUNCH'")

    payload: Dict[str, Any] = {
        "project_id": project_id,
        "starting_structure_id": starting_structure_id,
        "action": action,
    }

    if name:
        payload["name"] = name
    if description:
        payload["description"] = description
    if calculator_parameters:
        payload["calculator_parameters"] = calculator_parameters
    if calculator_kpoints:
        payload["calculator_kpoints"] = calculator_kpoints
    if calculator_basis_set:
        payload["calculator_basis_set"] = calculator_basis_set

    if extra_kwargs:
        payload.update(extra_kwargs)

    result = post(
        "api/vibes-simulation/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def get_vibes_simulation(simulation_id: str, **params: Any) -> Any:
    """
    Get a VIBES simulation by ID

    Args:
        simulation_id: str - The ID of the VIBES simulation
        **params: Additional GET parameters to pass to the API

    Returns:
        dict: The VIBES simulation data
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/vibes-simulation/{simulation_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    result = get(url)
    return result


def delete_vibes_simulation(simulation_id: str) -> Any:
    """
    Delete a VIBES simulation

    Args:
        simulation_id: str - The ID of the VIBES simulation to delete

    Returns:
        dict: The deletion response
    """
    result = delete(f"api/vibes-simulation/{simulation_id}/")
    return result


def associate_file_with_vibes_simulation(
    user_upload_id: str, simulation_id: str
) -> Any:
    """
    Associate a user upload file with a VIBES simulation

    Args:
        user_upload_id: str - The ID of the user upload
        simulation_id: str - The ID of the VIBES simulation

    Returns:
        dict: The file association data
    """
    result = post(
        "api/vibes-simulation-file/",
        payload={"user_upload": user_upload_id, "simulation": simulation_id},
    )
    return result


def get_vibes_simulation_files(simulation_id: str) -> Any:
    """
    Get files associated with a VIBES simulation

    Args:
        simulation_id: str - The ID of the VIBES simulation

    Returns:
        dict: List of files associated with the simulation
    """
    result = get(f"api/vibes-simulation-file/?simulation={simulation_id}")
    return result
