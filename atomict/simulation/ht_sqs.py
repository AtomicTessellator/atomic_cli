from typing import Any, Dict, List, Optional

from atomict.api import delete, get, post


def create_ht_sqs_exploration(
    project_id: str,
    name: str,
    target_concentrations: List[Dict[str, Any]],
    generated_permutations: List[Dict[str, float]],
    structure_id: str,
    num_structures_per_permutation: int = 1,
    structure_type: str = "userupload",
    description: Optional[str] = None,
    auto_max_size: bool = True,
    max_size: int = 8,
    atom_count_upper_limit: int = 200,
    cluster_cutoffs: Optional[List[float]] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Create a High Throughput SQS exploration

    Args:
        project_id: str - The project ID
        name: str - Name of the HT SQS exploration
        target_concentrations: list - List of {"element": str, "weight": float} dicts
        generated_permutations: list - List of concentration permutation dicts
        structure_id: str - ID of the structure to use as starting point
        num_structures_per_permutation: int - Number of SQS structures per permutation (default: 1)
        structure_type: str - Type of structure: "fhiaims", "mlrelax", or "userupload" (default: "userupload")
        description: str - Description of the exploration
        auto_max_size: bool - Auto-determine max supercell size (default: True)
        max_size: int - Maximum supercell size (default: 8)
        atom_count_upper_limit: int - Upper limit for atom count (default: 200)
        cluster_cutoffs: list - Cluster cutoff distances (default: [4.0, 4.0])
        extra_kwargs: dict - Additional parameters to pass to API

    Returns:
        dict: API response containing creation confirmation
    """

    # Validate structure_type
    valid_structure_types = ["fhiaims", "mlrelax", "userupload"]
    if structure_type not in valid_structure_types:
        raise ValueError(
            f"structure_type must be one of {valid_structure_types}, got '{structure_type}'"
        )

    # Validate target concentrations
    if not target_concentrations:
        raise ValueError("target_concentrations cannot be empty")

    for conc in target_concentrations:
        if not isinstance(conc, dict) or "element" not in conc or "weight" not in conc:
            raise ValueError(
                "Each target concentration must be a dict with 'element' and 'weight' keys"
            )
        if not 0 <= conc["weight"] <= 1.0:
            raise ValueError(
                f"Target concentration for {conc['element']} must be between 0 and 1"
            )

    # Check sum equals 1.0
    total_weight = sum(c["weight"] for c in target_concentrations)
    if abs(total_weight - 1.0) > 1e-6:
        raise ValueError(
            f"Sum of target concentrations must equal 1.0 (got {total_weight})"
        )

    # Validate generated_permutations
    if not generated_permutations:
        raise ValueError("generated_permutations cannot be empty")

    for i, perm in enumerate(generated_permutations):
        if not isinstance(perm, dict):
            raise ValueError(f"Permutation {i} must be a dict")

        # Check that each permutation's concentrations sum to 1.0
        perm_total = sum(perm.values())
        if abs(perm_total - 1.0) > 1e-6:
            raise ValueError(
                f"Permutation {i} concentrations must sum to 1.0 (got {perm_total})"
            )

    # Validate num_structures_per_permutation
    if num_structures_per_permutation < 1:
        raise ValueError("num_structures_per_permutation must be >= 1")

    # Set default cluster cutoffs
    if cluster_cutoffs is None:
        cluster_cutoffs = [4.0, 4.0]

    # Map structure_type to the correct API field
    structure_field_map = {
        "fhiaims": "starting_structure_id",
        "mlrelax": "starting_structure_mlrelax_id",
        "userupload": "starting_structure_userupload_id",
    }

    payload = {
        "project": project_id,
        "name": name,
        "description": description,
        "target_concentrations": target_concentrations,
        "generated_permutations": generated_permutations,
        "num_structures_per_permutation": num_structures_per_permutation,
        "auto_max_size": auto_max_size,
        "max_size": max_size,
        "atom_count_upper_limit": atom_count_upper_limit,
        "cluster_cutoffs": cluster_cutoffs,
    }

    # Add structure reference
    payload[structure_field_map[structure_type]] = structure_id

    if extra_kwargs:
        payload.update(extra_kwargs)

    result = post(
        "api/ht-sqs-exploration/",
        payload,
        extra_headers={"Content-Type": "application/json"},
    )

    return result


def get_ht_sqs_exploration(
    exploration_id: str, children: bool = False, **params: Any
) -> Dict[str, Any]:
    """
    Get a HT SQS exploration

    Args:
        exploration_id: str - The ID of the HT SQS exploration
        children: bool - Whether to include child SQS exploration details (default: False)
        **params: Additional GET parameters to pass to the API

    Returns:
        dict: HT SQS exploration details
    """
    # Start with any additional parameters
    query_params = params.copy()

    if children:
        query_params["children"] = "true"

    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    base_url = f"api/ht-sqs-exploration/{exploration_id}/"

    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url

    result = get(url)
    return result


def list_ht_sqs_explorations(
    project_id: Optional[str] = None,
    limit: Optional[int] = None,
    offset: Optional[int] = None,
    **params: Any,
) -> Dict[str, Any]:
    """
    List HT SQS explorations

    Args:
        project_id: str - Filter by project ID
        limit: int - Maximum number of results to return
        offset: int - Number of results to skip (for pagination)
        **params: Additional GET parameters to pass to the API

    Returns:
        dict: Paginated list of HT SQS explorations
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
    base_url = "api/ht-sqs-exploration/"

    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url

    result = get(url)
    return result


def delete_ht_sqs_exploration(exploration_id: str) -> Dict[str, Any]:
    """
    Delete a HT SQS exploration and all associated child SQS explorations

    Args:
        exploration_id: str - The ID of the HT SQS exploration to delete

    Returns:
        dict: API response confirming deletion
    """
    result = delete(f"api/ht-sqs-exploration/{exploration_id}/")
    return result


def get_ht_sqs_child(child_id: str, **params: Any) -> Dict[str, Any]:
    """
    Get a specific HT SQS child (HTSQS junction record)

    Args:
        child_id: str - The ID of the HTSQS child record
        **params: Additional GET parameters to pass to the API

    Returns:
        dict: HTSQS child details including linked SQS exploration
    """
    # Start with any additional parameters
    query_params = params.copy()

    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    base_url = f"api/ht-sqs/{child_id}/"

    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url

    result = get(url)
    return result


def list_ht_sqs_children(
    ht_sqs_exploration_id: Optional[str] = None,
    limit: Optional[int] = None,
    offset: Optional[int] = None,
    **params: Any,
) -> Dict[str, Any]:
    """
    List HT SQS children (HTSQS junction records)

    Args:
        ht_sqs_exploration_id: str - Filter by parent HT SQS exploration ID
        limit: int - Maximum number of results to return
        offset: int - Number of results to skip (for pagination)
        **params: Additional GET parameters to pass to the API

    Returns:
        dict: Paginated list of HTSQS children
    """
    query_params = params.copy()

    if ht_sqs_exploration_id:
        query_params["ht_sqs_exploration"] = ht_sqs_exploration_id
    if limit is not None:
        query_params["limit"] = limit
    if offset is not None:
        query_params["offset"] = offset

    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    base_url = "api/ht-sqs/"

    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url

    result = get(url)
    return result


# Constants for user convenience
DEFAULT_NUM_STRUCTURES_PER_PERMUTATION = 1
DEFAULT_STRUCTURE_TYPE = "userupload"
VALID_STRUCTURE_TYPES = ["fhiaims", "mlrelax", "userupload"]
DEFAULT_MAX_SIZE = 8
DEFAULT_ATOM_COUNT_UPPER_LIMIT = 200
DEFAULT_CLUSTER_CUTOFFS = [4.0, 4.0]
