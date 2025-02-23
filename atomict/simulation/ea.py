from atomict.api import get, post


def get_ea_exploration(exploration_id: str, **params):
    """
    Get EA exploration
    
    Args:
        exploration_id: str - The ID of the exploration
        **params: Additional GET parameters to pass to the API
    """
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration/{exploration_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def get_ea_exploration_sample(sample_id: str, **params):
    """
    Get EA exploration sample
    
    Args:
        sample_id: str - The ID of the sample
        **params: Additional GET parameters to pass to the API
    """
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration-sample/{sample_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def get_ea_exploration_samples(exploration_id: str, **params):
    """
    Get EA exploration samples
    
    Args:
        exploration_id: str - The ID of the exploration
        **params: Additional GET parameters to pass to the API
    """
    # Start with the required exploration parameter
    query_params = params.copy()
    query_params['exploration'] = exploration_id
    
    query_string = '&'.join(f"{k}={v}" for k, v in query_params.items())
    url = f"api/ea-exploration-sample/?{query_string}"
    return get(url)


def get_ea_exploration_analysis(analysis_id: str, **params):
    """
    Get EA exploration analysis
    
    Args:
        analysis_id: str - The ID of the analysis
        **params: Additional GET parameters to pass to the API
    """
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration-analysis/{analysis_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


def get_ea_exploration_analysis_file(analysis_file_id: str, **params):
    """
    Get EA exploration analysis file
    
    Args:
        analysis_file_id: str - The ID of the analysis file
        **params: Additional GET parameters to pass to the API
    """
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/ea-exploration-analysis-file/{analysis_file_id}/"
    url = f"{base_url}?{query_string}" if query_string else base_url
    return get(url)


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
