from typing import Optional
from atomict.api import delete, get, post


def get_simulation(simulation_id: str, full: bool = False, **params):
    """
    Get a SQS simulation

    Args:
        simulation_id: str - The ID of the simulation
        full: bool - Whether to get the full simulation details (default: False)
        **params: Additional GET parameters to pass to the API
    """
    # Start with any additional parameters
    query_params = params.copy()
    
    if full:
        query_params['full'] = 'true'
    
    # Build query string
    query_string = '&'.join(f"{k}={v}" for k, v in query_params.items())
    base_url = f"api/sqs-exploration/{simulation_id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def create_sqs_exploration(
    project_id: str,
    name: str,
    target_concentrations: list,
    structure_id: str,
    structure_type: str = "userupload",
    action: str = "DRAFT",
    description: Optional[str] = None,
    max_size: int = 8,
    cluster_cutoffs: Optional[list] = None,
    extra_exploration_kwargs: Optional[dict] = None,
):
    """
    Create a SQS exploration
    
    Args:
        project_id: str - The project ID
        name: str - Name of the exploration
        target_concentrations: list - List of {"element": str, "weight": float} dicts that must sum to 1.0
        structure_id: str - ID of the structure to use as starting point
        structure_type: str - Type of structure: "fhiaims", "mlrelax", or "userupload" (default: "userupload")
        action: str - "DRAFT" or "LAUNCH" (default: "DRAFT")
        description: str - Description of the exploration
        max_size: int - Maximum supercell size (default: 8)
        cluster_cutoffs: list - Cluster cutoff distances (default: [4.0, 4.0])
        extra_exploration_kwargs: dict - Additional parameters to pass to API
    
    Returns:
        dict: API response containing the created SQS exploration
    """
    
    if action not in ["DRAFT", "LAUNCH"]:
        raise ValueError("Action must be 'DRAFT' or 'LAUNCH'")
    
    # Validate structure_type
    valid_structure_types = ["fhiaims", "mlrelax", "userupload"]
    if structure_type not in valid_structure_types:
        raise ValueError(f"structure_type must be one of {valid_structure_types}, got '{structure_type}'")
    
    # Validate target concentrations
    if not target_concentrations:
        raise ValueError("target_concentrations cannot be empty")
    
    for conc in target_concentrations:
        if not isinstance(conc, dict) or "element" not in conc or "weight" not in conc:
            raise ValueError("Each target concentration must be a dict with 'element' and 'weight' keys")
        if not 0 <= conc["weight"] <= 1.0:
            raise ValueError(f"Target concentration for {conc['element']} must be between 0 and 1")
    
    # Check sum equals 1.0
    total_weight = sum(c["weight"] for c in target_concentrations)
    if abs(total_weight - 1.0) > 1e-6:
        raise ValueError(f"Sum of target concentrations must equal 1.0 (got {total_weight})")
    
    # Set default cluster cutoffs
    if cluster_cutoffs is None:
        cluster_cutoffs = [4.0, 4.0]
    
    payload = {
        "project": project_id,
        "name": name,
        "description": description,
        "target_concentrations": target_concentrations,
        "action": action,
        "max_size": max_size,
        "cluster_cutoffs": cluster_cutoffs,
    }
    
    # Map structure_type to the correct API field
    structure_field_map = {
        "fhiaims": "starting_structure_id",
        "mlrelax": "starting_structure_mlrelax_id", 
        "userupload": "starting_structure_userupload_id"
    }
    
    payload[structure_field_map[structure_type]] = structure_id
    
    if extra_exploration_kwargs:
        payload.update(extra_exploration_kwargs)
    
    result = post(
        "api/sqs-exploration/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )
    
    return result


def delete_sqs_exploration(exploration_id: str):
    """
    Delete a SQS exploration
    
    Args:
        exploration_id: str - The ID of the exploration to delete
        
    Returns:
        dict: API response confirming deletion
    """
    result = delete(f"api/sqs-exploration/{exploration_id}/")
    return result


def associate_user_upload_with_sqs_simulation(user_upload_id: str, exploration_id: str):
    """
    Associate a user upload with a SQS simulation
    """
    result = post(
        "api/sqs-simulation-file/",
        payload={"user_upload_id": user_upload_id, "exploration_id": exploration_id},
    )
    return result
