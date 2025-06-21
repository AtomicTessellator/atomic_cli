from typing import Dict, Any, Optional
from atomict.api import delete, get, patch, post
from atomict.simulation.models import (
    MODEL_ESEN_30M_OAM,
    MODEL_MATTERSIM_1_0_0_5M,
    MODEL_ORB_D3_V2,
    MODEL_ORB_V3_CONSERVATIVE,
)

# User-friendly model name mapping
ML_MODELS = {
    "orb_d3_v2": MODEL_ORB_D3_V2,
    "mattersim_1_0_0_5m": MODEL_MATTERSIM_1_0_0_5M,
    "orb_v3_conservative": MODEL_ORB_V3_CONSERVATIVE,
    "esen_30m_oam": MODEL_ESEN_30M_OAM,
}


def create_ht_mlrelax_exploration(
    project_id: str,
    name: str,
    source_ht_sqs_exploration_id: str,
    f_max: Optional[float] = None,
    model: str = "mattersim_1_0_0_5m",
    description: Optional[str] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Create a High Throughput ML Relaxation exploration

    This creates a batch of ML relaxation jobs from a completed HT SQS exploration.
    Individual ML relaxation tasks are automatically created for each completed 
    SQS structure in the source exploration.

    Args:
        project_id: str - The project ID
        name: str - Name of the HT ML relaxation exploration
        source_ht_sqs_exploration_id: str - ID of the completed HT SQS exploration
        f_max: float - Force convergence threshold (optional)
        model: str - ML model to use, one of:
            - "orb_d3_v2" (ORB D3 v2)
            - "mattersim_1_0_0_5m" (MatterSim 1.0.0.5M) [default]
            - "orb_v3_conservative" (ORB v3 Conservative OMAT)
            - "esen_30m_oam" (eSEN 30M OAM)
        description: str - Optional description
        extra_kwargs: dict - Additional parameters to pass to the API

    Returns:
        dict - API response containing creation status and count information

    Raises:
        ValueError: If model is not recognized
        APIValidationError: If source exploration is invalid or incomplete
    """
    # Validate model parameter
    if model not in ML_MODELS:
        valid_models = list(ML_MODELS.keys())
        raise ValueError(f"Invalid model '{model}'. Must be one of: {valid_models}")

    # Build payload with field mapping for LaunchableViewSet
    payload = {
        "project": project_id,  # LaunchableViewSet expects "project" not "project_id"
        "name": name,
        "source_ht_sqs_exploration": source_ht_sqs_exploration_id,
        "model": ML_MODELS[model],  # Convert user-friendly name to integer code
    }

    # Add optional parameters
    if f_max is not None:
        payload["f_max"] = f_max
    if description is not None:
        payload["description"] = description

    # Apply additional parameters
    if extra_kwargs:
        payload.update(extra_kwargs)

    # Call API endpoint
    result = post("api/ht-mlrelax/", payload)
    return result


def get_ht_mlrelax_exploration(exploration_id: str, include_children: bool = False, **params: Any) -> Dict[str, Any]:
    """
    Get High Throughput ML Relaxation exploration details

    Args:
        exploration_id: str - The HT ML relaxation exploration ID
        include_children: bool - Whether to include individual ML relaxation details
        **params: Additional GET parameters to pass to the API

    Returns:
        dict - HT ML relaxation exploration data with optional children
    """
    # Start with any additional parameters
    query_params = params.copy()
    
    if include_children:
        query_params["children"] = "true"
    
    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    base_url = f"api/ht-mlrelax/{exploration_id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def list_ht_mlrelax_explorations(
    project_id: Optional[str] = None,
    limit: Optional[int] = None,
    offset: Optional[int] = None,
    **params: Any,
) -> Dict[str, Any]:
    """
    List High Throughput ML Relaxation explorations

    Args:
        project_id: str - Filter by project ID
        limit: int - Maximum number of results to return
        offset: int - Number of results to skip (for pagination)
        **params: Additional GET parameters to pass to the API

    Returns:
        dict - Paginated list of HT ML relaxation explorations
    """
    query_params = params.copy()

    if project_id:
        query_params["project__id"] = project_id
    if limit is not None:
        query_params["limit"] = limit
    if offset is not None:
        query_params["offset"] = offset

    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    base_url = "api/ht-mlrelax/"

    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url

    result = get(url)
    return result


def update_ht_mlrelax_exploration(
    exploration_id: str,
    name: Optional[str] = None,
    f_max: Optional[float] = None,
    model: Optional[str] = None,
    description: Optional[str] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Update High Throughput ML Relaxation exploration parameters

    Args:
        exploration_id: str - The HT ML relaxation exploration ID
        name: str - New name for the exploration
        f_max: float - New force convergence threshold
        model: str - New ML model (same options as create_ht_mlrelax_exploration)
        description: str - New description
        extra_kwargs: dict - Additional parameters to pass to the API

    Returns:
        dict - Updated HT ML relaxation exploration data

    Raises:
        ValueError: If model is not recognized
    """
    payload = {}

    # Add optional parameters that are provided
    if name is not None:
        payload["name"] = name
    if f_max is not None:
        payload["f_max"] = f_max
    if description is not None:
        payload["description"] = description

    # Validate and convert model if provided
    if model is not None:
        if model not in ML_MODELS:
            valid_models = list(ML_MODELS.keys())
            raise ValueError(f"Invalid model '{model}'. Must be one of: {valid_models}")
        payload["model"] = ML_MODELS[model]

    # Apply additional parameters
    if extra_kwargs:
        payload.update(extra_kwargs)

    # Use PATCH for partial updates
    result = patch(f"api/ht-mlrelax/{exploration_id}/", payload)
    return result


def delete_ht_mlrelax_exploration(exploration_id: str) -> Dict[str, Any]:
    """
    Delete High Throughput ML Relaxation exploration and all associated child ML relaxations

    Args:
        exploration_id: str - The ID of the HT ML relaxation exploration to delete

    Returns:
        dict - API response confirming deletion
    """
    result = delete(f"api/ht-mlrelax/{exploration_id}/")
    return result


def launch_ht_mlrelax_exploration(
    exploration_id: str,
    cluster_id: Optional[str] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Launch High Throughput ML Relaxation exploration

    This transitions all individual ML relaxation tasks from DRAFT to READY status
    and submits them for execution on the specified cluster.

    Args:
        exploration_id: str - The HT ML relaxation exploration ID to launch
        cluster_id: str - Optional cluster ID for execution
        extra_kwargs: dict - Additional parameters to pass to the API

    Returns:
        dict - API response with launch confirmation

    Note:
        Launch functionality uses the LAUNCH action on the exploration.
        Individual ML relaxation tasks are automatically launched when the
        exploration is launched.
    """
    payload = {"action": "LAUNCH"}

    # Add cluster configuration if specified
    if cluster_id:
        if extra_kwargs is None:
            extra_kwargs = {}
        extra_kwargs["selected_cluster"] = cluster_id

    # Apply additional parameters
    if extra_kwargs:
        payload.update(extra_kwargs)

    # Use PATCH to trigger launch action
    result = patch(f"api/ht-mlrelax/{exploration_id}/", payload)
    return result
