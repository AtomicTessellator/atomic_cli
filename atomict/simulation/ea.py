from atomict.api import post
from atomict.resource_helpers import (
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_ea_exploration(exploration_id: str, **params):
    """Get EA exploration details.
    
    Args:
        exploration_id (str): The exploration identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/ea-exploration", exploration_id, **params)


def list_ea_explorations(**params):
    """List EA explorations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/ea-exploration", **params)


def create_ea_exploration(payload: dict[str, object]):
    """Create a new EA exploration.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/ea-exploration/", payload=payload)


def update_ea_exploration(exploration_id: str, fields: dict[str, object]):
    """Update EA exploration.
    
    Args:
        exploration_id (str): The exploration identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/ea-exploration", exploration_id, fields)


def delete_ea_exploration(exploration_id: str):
    """Delete EA exploration.
    
    Args:
        exploration_id (str): The exploration identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/ea-exploration", exploration_id)


def get_ea_exploration_sample(sample_id: str, **params):
    """Get EA exploration sample details.
    
    Args:
        sample_id (str): The sample identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/ea-exploration-sample", sample_id, **params)


def get_ea_exploration_samples(exploration_id: str, **params):
    """
    Get EA exploration samples
    
    Args:
        exploration_id: str - The ID of the exploration
        **params: Additional GET parameters to pass to the API
    """
    query_params = params.copy()
    query_params['exploration'] = exploration_id
    return list_resources("api/ea-exploration-sample", **query_params)


def list_ea_exploration_samples(**params):
    """List EA exploration samples.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/ea-exploration-sample", **params)


def update_ea_exploration_sample(sample_id: str, fields: dict[str, object]):
    """Update EA exploration sample.
    
    Args:
        sample_id (str): The sample identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/ea-exploration-sample", sample_id, fields)


def delete_ea_exploration_sample(sample_id: str):
    """Delete EA exploration sample.
    
    Args:
        sample_id (str): The sample identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/ea-exploration-sample", sample_id)


def get_ea_exploration_analysis(analysis_id: str, **params):
    """Get EA exploration analysis details.
    
    Args:
        analysis_id (str): The analysis identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/ea-exploration-analysis", analysis_id, **params)


def list_ea_exploration_analyses(**params):
    """List EA exploration analyses.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/ea-exploration-analysis", **params)


def create_ea_exploration_analysis(payload: dict[str, object]):
    """Create a new EA exploration analysis.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/ea-exploration-analysis/", payload=payload)


def update_ea_exploration_analysis(
    analysis_id: str, fields: dict[str, object]
):
    """Update EA exploration analysis.
    
    Args:
        analysis_id (str): The analysis identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/ea-exploration-analysis", analysis_id, fields)


def delete_ea_exploration_analysis(analysis_id: str):
    """Delete EA exploration analysis.
    
    Args:
        analysis_id (str): The analysis identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/ea-exploration-analysis", analysis_id)


def get_ea_exploration_analysis_file(analysis_file_id: str, **params):
    """Get EA exploration analysis file details.
    
    Args:
        analysis_file_id (str): The analysis file identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource(
        "api/ea-exploration-analysis-file", analysis_file_id, **params
    )


def list_ea_exploration_analysis_files(**params):
    """List EA exploration analysis files.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/ea-exploration-analysis-file", **params)


def create_ea_exploration_analysis_file(payload: dict[str, object]):
    """Create a new EA exploration analysis file.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/ea-exploration-analysis-file/", payload=payload)


def update_ea_exploration_analysis_file(
    file_id: str, fields: dict[str, object]
):
    """Update EA exploration analysis file.
    
    Args:
        file_id (str): The file identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/ea-exploration-analysis-file", file_id, fields)


def delete_ea_exploration_analysis_file(file_id: str):
    """Delete EA exploration analysis file.
    
    Args:
        file_id (str): The file identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/ea-exploration-analysis-file", file_id)


def associate_user_upload_with_ea_exploration(user_upload_id: str, analysis_id: str):
    return post(
        "api/ea-exploration-analysis-file/",
        payload={"user_upload_id": user_upload_id, "analysis_id": analysis_id},
    )


def create_exploration_sample(
    exploration_id: str,
    simulation_id: str = None,
    mlrelax_id: str = None,
    strain: float = None,
    matrix: int = None,
):
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


def create_ea_exploration_sample(
    exploration_id: str,
    simulation_id: str = None,
    mlrelax_id: str = None,
    strain: float = None,
    matrix: int = None,
):
    """Create a new EA exploration sample.
    
    Args:
        exploration_id (str): The exploration identifier.
        simulation_id (str | None): The simulation identifier.
        mlrelax_id (str | None): The ML relaxation identifier.
        strain (float | None): The strain value for the exploration sample.
        matrix (int | None): The matrix value for the exploration sample.
    
    Returns:
        dict: The API response payload.
    """
    return create_exploration_sample(
        exploration_id=exploration_id,
        simulation_id=simulation_id,
        mlrelax_id=mlrelax_id,
        strain=strain,
        matrix=matrix,
    )


def create_soec_exploration(payload: dict[str, object]):
    """Create a new SOEC exploration.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return create_ea_exploration(payload)


def get_soec_exploration(exploration_id: str, **params):
    """Get SOEC exploration details.
    
    Args:
        exploration_id (str): The exploration identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return get_ea_exploration(exploration_id, **params)


def list_soec_explorations(**params):
    """List SOEC explorations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_ea_explorations(**params)
