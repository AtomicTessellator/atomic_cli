from atomict.api import get, post


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


def associate_user_upload_with_sqs_simulation(user_upload_id: str, exploration_id: str):
    """
    Associate a user upload with a SQS simulation
    """
    result = post(
        "api/sqs-simulation-file/",
        payload={"user_upload_id": user_upload_id, "exploration_id": exploration_id},
    )
    return result
