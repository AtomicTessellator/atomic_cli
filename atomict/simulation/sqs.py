from atomict.api import post
from atomict.resource_helpers import (
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_simulation(simulation_id: str, full: bool = False, **params):
    """
    Get a SQS simulation

    Args:
        simulation_id: str - The ID of the simulation
        full: bool - Whether to get the full simulation details (default: False)
        **params: Additional GET parameters to pass to the API
    """
    query_params = params.copy()
    if full:
        query_params["full"] = "true"
    return retrieve_resource("api/sqs-exploration", simulation_id, **query_params)


def get_sqs_exploration(exploration_id: str, **params):
    """Get SQS exploration details.
    
    Args:
        exploration_id (str): The exploration identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return get_simulation(exploration_id, **params)


def list_sqs_explorations(**params):
    """List SQS explorations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/sqs-exploration", **params)


def create_sqs_exploration(payload: dict[str, object]):
    """Create a new SQS exploration.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/sqs-exploration/", payload=payload)


def update_sqs_exploration(exploration_id: str, fields: dict[str, object]):
    """Update SQS exploration.
    
    Args:
        exploration_id (str): The exploration identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/sqs-exploration", exploration_id, fields)


def delete_sqs_exploration(exploration_id: str):
    """Delete SQS exploration.
    
    Args:
        exploration_id (str): The exploration identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/sqs-exploration", exploration_id)


def associate_user_upload_with_sqs_simulation(user_upload_id: str, exploration_id: str):
    """
    Associate a user upload with a SQS simulation
    """
    result = post(
        "api/sqs-simulation-file/",
        payload={"user_upload_id": user_upload_id, "exploration_id": exploration_id},
    )
    return result


def get_sqs_target_concentration(target_concentration_id: str, **params):
    """Get SQS target concentration entry details.
    
    Args:
        target_concentration_id (str): The target concentration identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource(
        "api/sqs-target-concentration", target_concentration_id, **params
    )


def list_sqs_target_concentrations(**params):
    """List target concentrations for an SQS exploration.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/sqs-target-concentration", **params)


def create_sqs_target_concentration(payload: dict[str, object]):
    """Create a new SQS target concentration entry.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/sqs-target-concentration/", payload=payload)


def update_sqs_target_concentration(
    target_concentration_id: str, fields: dict[str, object]
):
    """Update SQS target concentration entry.
    
    Args:
        target_concentration_id (str): The target concentration identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource(
        "api/sqs-target-concentration", target_concentration_id, fields
    )


def delete_sqs_target_concentration(target_concentration_id: str):
    """Delete SQS target concentration entry.
    
    Args:
        target_concentration_id (str): The target concentration identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/sqs-target-concentration", target_concentration_id)


def get_sqs_simulation_file(file_id: str, **params):
    """Get SQS simulation file details.
    
    Args:
        file_id (str): The file identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/sqs-simulation-file", file_id, **params)


def list_sqs_simulation_files(**params):
    """List SQS simulation files for an exploration.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/sqs-simulation-file", **params)


def create_sqs_simulation_file(payload: dict[str, object]):
    """Create a new SQS simulation file.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/sqs-simulation-file/", payload=payload)


def update_sqs_simulation_file(file_id: str, fields: dict[str, object]):
    """Update SQS simulation file.
    
    Args:
        file_id (str): The file identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/sqs-simulation-file", file_id, fields)


def delete_sqs_simulation_file(file_id: str):
    """Delete SQS simulation file.
    
    Args:
        file_id (str): The file identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/sqs-simulation-file", file_id)
