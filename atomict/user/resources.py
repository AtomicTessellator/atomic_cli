from atomict.api import get
from atomict.resource_helpers import delete_resource, list_resources, retrieve_resource, update_resource


def get_user(user_id: str, **params) -> dict:
    """Get user details.
    
    Args:
        user_id (str): The user identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/user", user_id, **params)


def list_users(**params) -> dict:
    """List users.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/user", **params)


def update_user(user_id: str, fields: dict[str, object]) -> dict:
    """Update user.
    
    Args:
        user_id (str): The user identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/user", user_id, fields)


def delete_user(user_id: str) -> dict:
    """Delete user.
    
    Args:
        user_id (str): The user identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/user", user_id)


def lookup_simulation_workspace(simulation_uuid: str, **params) -> dict:
    """Look up the workspace for a simulation.
    
    Args:
        simulation_uuid (str): The simulation UUID.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    query = dict(params)
    return get(f"simulation/lookup/{simulation_uuid}/", params=query or None)
