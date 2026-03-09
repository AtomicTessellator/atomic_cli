from atomict.api import get, post
from atomict.simulation.models import MODEL_ORB_D3_V2, MODEL_MATTERSIM_1_0_0_5M, MODEL_ORB_V3_CONSERVATIVE, MODEL_ESEN_30M_OAM


def get_phonon_run(id: str, **params):
    """Get phonon run details.
    
    Args:
        id (str): The resource identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    # Build query string from parameters
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/phonon-run/{id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def get_phonon_sim_run(id: str, **params):
    """Get phonon simulation run details.
    
    Args:
        id (str): The resource identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    # Build query string from parameters
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/phonon-run-simulation/{id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def associate_user_upload_with_phonon_sim_run(user_upload_id: str, phonon_run_id: str):
    """Associate a user upload with a phonon simulation run.
    
    Args:
        user_upload_id (str): The user upload identifier.
        phonon_run_id (str): The phonon run identifier.
    
    Returns:
        dict: The API response payload.
    """
    result = post(
        "api/phonon-run-simulation-file/",
        payload={"user_upload_id": user_upload_id, "phono3py_run_simulation_id": phonon_run_id},
    )
    return result


def create_phonon_run(project_id: str, source_geometry_id: str, action: str, name: str = None, description: str = None, model: int = MODEL_ORB_D3_V2, extra_simulation_kwargs: dict = None):
    """Create a new phonon run.
    
    Args:
        project_id (str): The project identifier.
        source_geometry_id (str): The source geometry identifier.
        action (str): The action to perform for the simulation.
        name (str | None): The resource name.
        description (str | None): The resource description.
        model (int): The model identifier to use for the run.
        extra_simulation_kwargs (dict | None): Additional simulation keyword arguments to include in the request.
    
    Returns:
        dict: The API response payload.
    """

    if action not in ["LAUNCH", "DRAFT"]:
        raise ValueError(f"Invalid action: {action} (must be 'LAUNCH' or 'DRAFT')")
    
    # Validate model is one of the supported constants
    valid_models = [
        MODEL_ORB_D3_V2,
        MODEL_MATTERSIM_1_0_0_5M,
        MODEL_ORB_V3_CONSERVATIVE,
        MODEL_ESEN_30M_OAM,
    ]
    if model not in valid_models:
        raise ValueError(f"Invalid model: {model}")
    
    payload = {
        "project_id": project_id,
        "source_geometry_id": source_geometry_id,
        "action": action,
        "name": name,
        "description": description,
        "model": model,
        "extra_simulation_kwargs": extra_simulation_kwargs,
    }
    result = post("api/phonon-run/", payload=payload)
    return result


def get_phonon_sim_run_files(phonon_sim_run_id: str):
    """List phonon simulation run files.
    
    Args:
        phonon_sim_run_id (str): The phonon simulation run identifier.
    
    Returns:
        dict: The API response payload.
    """
    result = get(f"api/phonon-run-simulation-file/?phono3py_run_simulation__id={phonon_sim_run_id}")
    return result
