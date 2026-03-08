from atomict.api import post
from atomict.resource_helpers import (
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_ea_exploration(exploration_id: str, **params):
    """
    Get EA exploration
    
    Args:
        exploration_id: str - The ID of the exploration
        **params: Additional GET parameters to pass to the API
    """
    return retrieve_resource("api/ea-exploration", exploration_id, **params)


def list_ea_explorations(**params):
    return list_resources("api/ea-exploration", **params)


def create_ea_exploration(payload: dict[str, object]):
    return post("api/ea-exploration/", payload=payload)


def update_ea_exploration(exploration_id: str, fields: dict[str, object]):
    return update_resource("api/ea-exploration", exploration_id, fields)


def delete_ea_exploration(exploration_id: str):
    return delete_resource("api/ea-exploration", exploration_id)


def get_ea_exploration_sample(sample_id: str, **params):
    """
    Get EA exploration sample
    
    Args:
        sample_id: str - The ID of the sample
        **params: Additional GET parameters to pass to the API
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
    return list_resources("api/ea-exploration-sample", **params)


def update_ea_exploration_sample(sample_id: str, fields: dict[str, object]):
    return update_resource("api/ea-exploration-sample", sample_id, fields)


def delete_ea_exploration_sample(sample_id: str):
    return delete_resource("api/ea-exploration-sample", sample_id)


def get_ea_exploration_analysis(analysis_id: str, **params):
    """
    Get EA exploration analysis
    
    Args:
        analysis_id: str - The ID of the analysis
        **params: Additional GET parameters to pass to the API
    """
    return retrieve_resource("api/ea-exploration-analysis", analysis_id, **params)


def list_ea_exploration_analyses(**params):
    return list_resources("api/ea-exploration-analysis", **params)


def create_ea_exploration_analysis(payload: dict[str, object]):
    return post("api/ea-exploration-analysis/", payload=payload)


def update_ea_exploration_analysis(
    analysis_id: str, fields: dict[str, object]
):
    return update_resource("api/ea-exploration-analysis", analysis_id, fields)


def delete_ea_exploration_analysis(analysis_id: str):
    return delete_resource("api/ea-exploration-analysis", analysis_id)


def get_ea_exploration_analysis_file(analysis_file_id: str, **params):
    """
    Get EA exploration analysis file
    
    Args:
        analysis_file_id: str - The ID of the analysis file
        **params: Additional GET parameters to pass to the API
    """
    return retrieve_resource(
        "api/ea-exploration-analysis-file", analysis_file_id, **params
    )


def list_ea_exploration_analysis_files(**params):
    return list_resources("api/ea-exploration-analysis-file", **params)


def create_ea_exploration_analysis_file(payload: dict[str, object]):
    return post("api/ea-exploration-analysis-file/", payload=payload)


def update_ea_exploration_analysis_file(
    file_id: str, fields: dict[str, object]
):
    return update_resource("api/ea-exploration-analysis-file", file_id, fields)


def delete_ea_exploration_analysis_file(file_id: str):
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
    return create_exploration_sample(
        exploration_id=exploration_id,
        simulation_id=simulation_id,
        mlrelax_id=mlrelax_id,
        strain=strain,
        matrix=matrix,
    )


def create_soec_exploration(payload: dict[str, object]):
    return create_ea_exploration(payload)


def get_soec_exploration(exploration_id: str, **params):
    return get_ea_exploration(exploration_id, **params)


def list_soec_explorations(**params):
    return list_ea_explorations(**params)
