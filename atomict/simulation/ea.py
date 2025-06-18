from typing import Any, Dict, List, Optional
from atomict.api import delete, get, post


def get_ea_exploration(exploration_id: str, **params) -> Dict[str, Any]:
    """
    Get EA exploration

    Args:
        exploration_id: str - The ID of the exploration
        **params: Additional GET parameters to pass to the API
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration/{exploration_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def get_ea_exploration_sample(sample_id: str, **params) -> Dict[str, Any]:
    """
    Get EA exploration sample

    Args:
        sample_id: str - The ID of the sample
        **params: Additional GET parameters to pass to the API
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration-sample/{sample_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def get_ea_exploration_samples(exploration_id: str, **params) -> Dict[str, Any]:
    """
    Get EA exploration samples

    Args:
        exploration_id: str - The ID of the exploration
        **params: Additional GET parameters to pass to the API
    """
    # Start with the required exploration parameter
    query_params = params.copy()
    query_params["exploration"] = exploration_id

    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    url = f"api/ea-exploration-sample/?{query_string}"
    return get(url)


def get_ea_exploration_analysis(analysis_id: str, **params) -> Dict[str, Any]:
    """
    Get EA exploration analysis

    Args:
        analysis_id: str - The ID of the analysis
        **params: Additional GET parameters to pass to the API
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration-analysis/{analysis_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def get_ea_exploration_analysis_file(analysis_file_id: str, **params) -> Dict[str, Any]:
    """
    Get EA exploration analysis file

    Args:
        analysis_file_id: str - The ID of the analysis file
        **params: Additional GET parameters to pass to the API
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration-analysis-file/{analysis_file_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def associate_user_upload_with_ea_exploration(
    user_upload_id: str, analysis_id: str
) -> Dict[str, Any]:
    return post(
        "api/ea-exploration-analysis-file/",
        payload={"user_upload_id": user_upload_id, "analysis_id": analysis_id},
    )


# Constants for EA exploration configuration
STRESS_ALGORITHMS = {"ULICS": 0, "OHESS": 1, "ASESS": 2}

STRESS_METHODS = {"Static": 0, "Dynamic": 1}

CALCULATORS = {
    "FHI-AIMS": 0,
    "MATTER_SIM": 1,  # Deprecated
    "ORB_D3": 2,  # Deprecated
    "ORB_V3_CONSERVATIVE": 3,
    "ESEN_30M_OAM": 4,
}


def create_ea_exploration(
    project_id: str,
    name: str,
    structure_id: str,
    structure_type: str = "userupload",
    action: str = "DRAFT",
    strains_list: Optional[List[float]] = None,
    stress_algorithm: int = 2,  # ASESS default
    stress_method: int = 1,  # Dynamic default
    calculator: int = 0,  # FHI-AIMS default
    num_last_samples: int = 1000,
    description: Optional[str] = None,
    make_conventional_cell: bool = False,
    remove_spurious_distortions: bool = True,
    add_vacuum: Optional[int] = None,
    is_ht: bool = False,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Create an EA (SOEC) exploration for elastic analysis

    Args:
        project_id: str - The project ID
        name: str - Name of the exploration
        structure_id: str - ID of the structure to use as starting point
        structure_type: str - Type of structure: "fhiaims", "mlrelax", or
            "userupload" (default: "userupload")
        action: str - "DRAFT" or "LAUNCH" (default: "DRAFT")
        strains_list: list - List of strain values
            (default: [-0.06, -0.03, 0.03, 0.06])
        stress_algorithm: int - Stress algorithm: 0=ULICS, 1=OHESS, 2=ASESS
            (default: 2)
        stress_method: int - Stress method: 0=Static, 1=Dynamic (default: 1)
        calculator: int - Calculator type: 0=FHI-AIMS, 3=ORB_V3, 4=eSEN
            (default: 0)
        num_last_samples: int - Number of last samples to use (default: 1000)
        description: str - Description of the exploration
        make_conventional_cell: bool - Whether to make conventional cell
            (default: False)
        remove_spurious_distortions: bool - Whether to remove spurious
            distortions (default: True)
        add_vacuum: int - Amount of vacuum to add in Å (0-100, default: None)
        is_ht: bool - High-throughput mode (default: False)
        extra_kwargs: dict - Additional parameters to pass to API

    Returns:
        dict: API response containing the created EA exploration
    """

    if action not in ["DRAFT", "LAUNCH"]:
        raise ValueError("Action must be 'DRAFT' or 'LAUNCH'")

    # Validate structure_type
    valid_structure_types = ["fhiaims", "mlrelax", "userupload"]
    if structure_type not in valid_structure_types:
        raise ValueError(
            f"structure_type must be one of {valid_structure_types}, "
            f"got '{structure_type}'"
        )

    # Set default strains list
    if strains_list is None:
        strains_list = [-0.06, -0.03, 0.03, 0.06]

    # Validate strains list
    if len(strains_list) < 4:
        raise ValueError("strains_list must contain at least 4 elements")

    # Validate add_vacuum range
    if add_vacuum is not None and (add_vacuum < 0 or add_vacuum > 100):
        raise ValueError("add_vacuum must be between 0 and 100 Å")

    payload = {
        "project": project_id,  # LaunchableViewSet expects "project" not "project_id"
        "name": name,
        "description": description,
        "action": action,
        "strains_list": strains_list,
        "stress_algorithm": stress_algorithm,
        "stress_method": stress_method,
        "calculator": calculator,
        "num_last_samples": num_last_samples,
        "make_conventional_cell": make_conventional_cell,
        "remove_spurious_distortions": remove_spurious_distortions,
        "add_vacuum": add_vacuum,
        "is_ht": is_ht,
    }

    # Map structure_type to the correct API field
    structure_field_map = {
        "fhiaims": "starting_structure_id",
        "mlrelax": "starting_structure_mlrelax_id",
        "userupload": "starting_structure_userupload_id",
    }

    payload[structure_field_map[structure_type]] = structure_id

    if extra_kwargs:
        payload.update(extra_kwargs)

    result = post(
        "api/ea-exploration/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def delete_ea_exploration(exploration_id: str) -> Dict[str, Any]:
    """
    Delete an EA exploration

    Args:
        exploration_id: str - The ID of the exploration to delete

    Returns:
        dict: API response confirming deletion
    """
    result = delete(f"api/ea-exploration/{exploration_id}/")
    return result


def delete_ea_exploration_sample(sample_id: str) -> Dict[str, Any]:
    """
    Delete an EA exploration sample

    Args:
        sample_id: str - The ID of the sample to delete

    Returns:
        dict: API response confirming deletion
    """
    result = delete(f"api/ea-exploration-sample/{sample_id}/")
    return result


def create_ea_exploration_analysis(
    exploration_id: str,
    compute_directional_properties: bool = True,
    action: str = "DRAFT",
) -> Dict[str, Any]:
    """
    Create an EA exploration analysis

    Args:
        exploration_id: str - The ID of the exploration
        compute_directional_properties: bool - Whether to compute directional
            properties (default: True)
        action: str - "DRAFT" or "LAUNCH" (default: "DRAFT")

    Returns:
        dict: API response containing the created analysis
    """

    if action not in ["DRAFT", "LAUNCH"]:
        raise ValueError("Action must be 'DRAFT' or 'LAUNCH'")

    payload = {
        "exploration_id": exploration_id,
        "compute_directional_properties": compute_directional_properties,
        "action": action,
    }

    result = post(
        "api/ea-exploration-analysis/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def delete_ea_exploration_analysis(analysis_id: str) -> Dict[str, Any]:
    """
    Delete an EA exploration analysis

    Args:
        analysis_id: str - The ID of the analysis to delete

    Returns:
        dict: API response confirming deletion
    """
    result = delete(f"api/ea-exploration-analysis/{analysis_id}/")
    return result


def create_ea_exploration_analysis_file(
    analysis_id: str, user_upload_id: str
) -> Dict[str, Any]:
    """
    Create an EA exploration analysis file by associating a user upload

    Args:
        analysis_id: str - The ID of the analysis
        user_upload_id: str - The ID of the user upload to associate

    Returns:
        dict: API response containing the created analysis file
    """
    payload = {
        "analysis_id": analysis_id,
        "user_upload_id": user_upload_id,
    }

    result = post(
        "api/ea-exploration-analysis-file/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def delete_ea_exploration_analysis_file(analysis_file_id: str) -> Dict[str, Any]:
    """
    Delete an EA exploration analysis file

    Args:
        analysis_file_id: str - The ID of the analysis file to delete

    Returns:
        dict: API response confirming deletion
    """
    result = delete(f"api/ea-exploration-analysis-file/{analysis_file_id}/")
    return result


def create_exploration_sample(
    exploration_id: str,
    simulation_id: Optional[str] = None,
    mlrelax_id: Optional[str] = None,
    strain: Optional[float] = None,
    matrix: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Create an exploration sample

    exploration_id: str - The ID of the exploration to associate the sample with
    simulation_id: str - The ID of the simulation to associate with the exploration
    strain: float - The strain to associate with the sample
    matrix: int - The matrix to associate with the sample
    """

    if simulation_id is None and mlrelax_id is None:
        raise ValueError("Either simulation_id or mlrelax_id must be provided")

    payload = {
        "exploration_id": exploration_id,
        "strain": strain,
        "matrix": matrix,
    }

    if simulation_id:
        payload["simulation_id"] = simulation_id
    elif mlrelax_id:
        payload["mlrelax_id"] = mlrelax_id

    return post(
        "api/ea-exploration-sample/",
        payload=payload,
    )
