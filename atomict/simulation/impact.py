from atomict.api import get, post


def get_impact_simulation(id: str, **params):
    """Get impact simulation details.
    
    Args:
        id (str): The resource identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    # Build query string from parameters
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/impact-simulation/{id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def associate_user_upload_with_impact_simulation(user_upload_id: str, impact_simulation_id: str):
    """Associate a user upload with an impact simulation.
    
    Args:
        user_upload_id (str): The user upload identifier.
        impact_simulation_id (str): The impact simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    result = post(
        "api/impact-simulation-file/",
        payload={"user_upload_id": user_upload_id, "impact_simulation_id": impact_simulation_id},
    )
    return result


def get_impact_simulation_files(impact_simulation_id: str):
    """List files associated with an impact simulation.
    
    Args:
        impact_simulation_id (str): The impact simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    result = get(f"api/impact-simulation-file/?impact_simulation__id={impact_simulation_id}")
    return result
