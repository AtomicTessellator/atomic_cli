from atomict.api import get, post, patch, delete
from atomict.exceptions import APIValidationError


# Structure type to API field mapping
STRUCTURE_FIELD_MAP = {
    "fhiaims": "starting_structure_id",
    "mlrelax": "starting_structure_mlrelax_id",
    "userupload": "starting_structure_userupload_id",
}


def create_kpoint_exploration(
    project_id: str,
    name: str,
    structure_id: str,
    structure_type: str = "userupload",
    k_point_range_lower: int = 3,
    k_point_range_upper: int = 6,
    evenly_spaced_kpoints: bool = False,
    action: str = "DRAFT",
    extra_kwargs: dict | None = None,
):
    """
    Create a K-point convergence exploration

    Args:
        project_id: Project ID to create exploration in
        name: Name for the exploration
        structure_id: ID of the starting structure
        structure_type: Type of structure ("fhiaims", "mlrelax", "userupload")
        k_point_range_lower: Lower bound of k-point range (default: 3)
        k_point_range_upper: Upper bound of k-point range (default: 6)
        evenly_spaced_kpoints: Whether to use evenly spaced k-points (default: False)
        action: "DRAFT" or "LAUNCH"
        extra_kwargs: Additional parameters for the API call

    Returns:
        API response dictionary
    """
    # Validate action
    if action not in ["DRAFT", "LAUNCH"]:
        raise APIValidationError(
            f"Invalid action '{action}'. Must be 'DRAFT' or 'LAUNCH'"
        )

    # Validate structure type
    if structure_type not in STRUCTURE_FIELD_MAP:
        valid_types = list(STRUCTURE_FIELD_MAP.keys())
        raise APIValidationError(
            f"Invalid structure_type '{structure_type}'. Must be one of: {valid_types}"
        )

    # Validate k-point range
    if k_point_range_lower >= k_point_range_upper:
        raise APIValidationError(
            "k_point_range_lower must be less than k_point_range_upper"
        )

    # Build payload with field mapping
    payload = {
        "project": project_id,
        "name": name,
        "k_point_range_lower": k_point_range_lower,
        "k_point_range_upper": k_point_range_upper,
        "evenly_spaced_kpoints": evenly_spaced_kpoints,
        "action": action,
        STRUCTURE_FIELD_MAP[structure_type]: structure_id,
    }

    # Add extra parameters if provided
    if extra_kwargs:
        payload.update(extra_kwargs)

    result = post("api/kpoint-exploration/", payload=payload)
    return result


def delete_kpoint_exploration(exploration_id: str):
    """
    Delete a K-point exploration

    Args:
        exploration_id: ID of the exploration to delete

    Returns:
        API response dictionary
    """
    result = delete(f"api/kpoint-exploration/{exploration_id}/")
    return result


def get_kpoint_exploration(simulation_id: str):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-exploration/{simulation_id}/")
    return result


def update_kpoint_exploration(exploration_id: str, fields: dict):
    """
    Update KPoint exploration
    """
    result = patch(f"api/kpoint-exploration/{exploration_id}/", payload=fields)
    return result


def get_kpoint_simulation_list(exploration_id: str):
    """
    Get kpoints for a simulation
    """
    result = get(f"api/kpoint-simulation/?exploration__id={exploration_id}")
    return result


def get_kpoint_simulation(simulation_id: str):
    """
    Get KPoint simulation
    """
    result = get(f"api/kpoint-simulation/{simulation_id}/")
    return result


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


def get_kpoint_analysis(analysis_id: str):
    """
    Get KPoint analysis
    """
    result = get(f"api/kpoint-analysis/{analysis_id}/")
    return result


def update_kpoint_analysis(analysis_id: str, fields: dict):
    """
    Update KPoint analysis
    """
    result = patch(f"api/kpoint-analysis/{analysis_id}/", payload=fields)
    return result


def delete_kpoint_simulation(simulation_id: str):
    """
    Delete a K-point simulation

    Args:
        simulation_id: ID of the simulation to delete
    """
    result = delete(f"api/kpoint-simulation/{simulation_id}/")
    return result


def create_kpoint_analysis(exploration_id: str, extra_kwargs: dict | None = None):
    """
    Create a K-point analysis for an exploration

    Args:
        exploration_id: ID of the exploration to analyze
        extra_kwargs: Additional parameters for the API call
    """
    payload = {"exploration_id": exploration_id}

    # Add extra parameters if provided
    if extra_kwargs:
        payload.update(extra_kwargs)

    result = post("api/kpoint-analysis/", payload=payload)
    return result


def delete_kpoint_analysis(analysis_id: str):
    """
    Delete a K-point analysis

    Args:
        analysis_id: ID of the analysis to delete
    """
    result = delete(f"api/kpoint-analysis/{analysis_id}/")
    return result
