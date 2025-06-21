from typing import Any, Dict, Optional

from atomict.api import delete, get, patch, post

# Example queries for user guidance
EXAMPLE_QUERIES = {
    "formula": "formula=Cu2O",
    "elements": "Cu,O",
    "space_group": "spacegroup=225",
    "energy_range": "energy<-3.0",
    "combined": "formula=Cu2O,spacegroup=225",
    "structure_count": "natoms>10,natoms<50",
    "magnetic": "magmom_per_atom>0.1",
}


def validate_element_symbol(element: str) -> bool:
    """
    Validate chemical element symbol

    Args:
        element: Chemical element symbol to validate

    Returns:
        bool: True if valid element symbol
    """
    # Basic validation - should be 1-2 characters, first uppercase
    if not element or len(element) > 2:
        return False

    if not element[0].isupper():
        return False

    if len(element) == 2 and not element[1].islower():
        return False

    return True


def validate_ase_query(query: str) -> bool:
    """
    Basic ASE query syntax validation

    Args:
        query: ASE database query string

    Returns:
        bool: True if query appears valid
    """
    if not query or not isinstance(query, str):
        return False

    # Basic checks for common query patterns
    return True  # For now, let the backend handle detailed validation


def search_structures(
    query: str,
    limit: int = 30,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Search structures using ASE database query syntax

    Requires utils dependencies: uv pip install -e ".[utils]"

    Args:
        query: ASE database query string (e.g., "formula=Cu2O")
        limit: Maximum results (1-30, default: 30)
        extra_kwargs: Additional GET parameters

    Returns:
        dict: {"status": "OK", "results": [structure_dicts...]}

    Example:
        >>> search_structures("formula=Cu2O")
        >>> search_structures("spacegroup=225", limit=10)
    """
    if not validate_ase_query(query):
        raise ValueError("Invalid ASE query format")

    if not isinstance(limit, int) or limit < 1 or limit > 30:
        raise ValueError("Limit must be an integer between 1 and 30")

    # Build query parameters
    params = {"query": query}

    if extra_kwargs:
        params.update(extra_kwargs)

    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"structure-search/?{query_string}"

    return get(url)


def search_structures_by_element(
    element: str,
    exp_observed: bool = True,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Search structures containing specific element from Materials Project data

    Args:
        element: Chemical element symbol (e.g., "Cu", "Fe", "O")
        exp_observed: Filter for experimentally observed structures (default: True)
        extra_kwargs: Additional GET parameters

    Returns:
        dict: {"results": [{"mp_id", "cif_structure", "ase_atoms", "data", ...}]}

    Example:
        >>> search_structures_by_element("Cu")
        >>> search_structures_by_element("Fe", exp_observed=False)
    """
    if not validate_element_symbol(element):
        raise ValueError(f"Invalid element symbol: '{element}'")

    # Build query parameters
    params = {"element": element}

    if not exp_observed:
        params["exp_observed"] = "false"

    if extra_kwargs:
        params.update(extra_kwargs)

    # Build query string
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"reality/matprj/element-structures/?{query_string}"

    return get(url)


def create_discovery_query(
    project_id: str,
    name: str,
    query: str,
    description: Optional[str] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Create discovery query for database exploration

    Args:
        project_id: The project ID to associate with
        name: Name for the discovery query
        query: SQL query string for database exploration
        description: Optional description of the query
        extra_kwargs: Additional parameters

    Returns:
        dict: Created discovery query object

    Example:
        >>> create_discovery_query("proj123", "Find Cu compounds",
        ...                       "SELECT * FROM structures WHERE formula LIKE '%Cu%'")
    """
    if not project_id:
        raise ValueError("project_id is required")
    if not name:
        raise ValueError("name is required")
    if not query:
        raise ValueError("query is required")

    payload = {
        "project": project_id,  # Use "project" not "project_id" for consistency
        "name": name,
        "query": query,
        "description": description,
    }

    if extra_kwargs:
        payload.update(extra_kwargs)

    return post("api/discovery-query/", payload=payload)


def get_discovery_query(query_id: str, **params) -> Dict[str, Any]:
    """
    Get discovery query by ID

    Args:
        query_id: The ID of the discovery query
        **params: Additional GET parameters

    Returns:
        dict: Discovery query object
    """
    query_string = "&".join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/discovery-query/{query_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url

    return get(url)


def list_discovery_queries(
    project_id: Optional[str] = None, **params
) -> Dict[str, Any]:
    """
    List discovery queries, optionally filtered by project

    Args:
        project_id: Optional project ID to filter by
        **params: Additional GET parameters

    Returns:
        dict: List of discovery query objects

    Example:
        >>> list_discovery_queries(project_id="proj123")
        >>> list_discovery_queries()  # All queries user has access to
    """
    query_params = params.copy()

    if project_id:
        query_params["project__id"] = project_id

    query_string = "&".join(f"{k}={v}" for k, v in query_params.items())
    url = (
        f"api/discovery-query/?{query_string}"
        if query_string
        else "api/discovery-query/"
    )

    return get(url)


def update_discovery_query(
    query_id: str,
    name: Optional[str] = None,
    query: Optional[str] = None,
    description: Optional[str] = None,
    extra_kwargs: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Update discovery query

    Args:
        query_id: The ID of the discovery query to update
        name: Optional new name
        query: Optional new query string
        description: Optional new description
        extra_kwargs: Additional parameters

    Returns:
        dict: Updated discovery query object
    """
    payload = {}

    if name is not None:
        payload["name"] = name
    if query is not None:
        payload["query"] = query
    if description is not None:
        payload["description"] = description

    if extra_kwargs:
        payload.update(extra_kwargs)

    if not payload:
        raise ValueError("At least one field must be provided for update")

    return patch(f"api/discovery-query/{query_id}/", payload=payload)


def delete_discovery_query(query_id: str) -> Dict[str, Any]:
    """
    Delete discovery query

    Args:
        query_id: The ID of the discovery query to delete

    Returns:
        dict: Deletion confirmation
    """
    return delete(f"api/discovery-query/{query_id}/")
