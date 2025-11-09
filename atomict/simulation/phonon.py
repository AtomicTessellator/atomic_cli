from atomict.api import get, post
from atomict.simulation.models import MODEL_ORB_D3_V2, MODEL_MATTERSIM_1_0_0_5M, MODEL_ORB_V3_CONSERVATIVE, MODEL_ESEN_30M_OAM


def get_phonon_run(id: str, **params):
    """
    Get Phonon Run
    
    A Phonon Run is a collection of Phonon Simulations.

    Args:
        id: str - The ID of the Phonon Run
        **params: Additional GET parameters to pass to the API
    """
    # Build query string from parameters
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/phonon-run/{id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def get_phonon_sim_run(id: str, **params):
    """
    Get Phonon Simulation Run

    Args:
        id: str - The ID of the Phonon Simulation Run
        **params: Additional GET parameters to pass to the API
    """
    # Build query string from parameters
    query_string = '&'.join(f"{k}={v}" for k, v in params.items())
    base_url = f"api/phonon-run-simulation/{id}/"
    
    # Add query string if we have parameters
    url = f"{base_url}?{query_string}" if query_string else base_url
    
    result = get(url)
    return result


def associate_user_upload_with_phonon_sim_run(user_upload_id: str, phonon_run_id: str):
    """
    Associate a user upload with a Phonon Simulation Run
    """
    result = post(
        "api/phonon-run-simulation-file/",
        payload={"user_upload_id": user_upload_id, "phono3py_run_simulation_id": phonon_run_id},
    )
    return result


def create_phonon_run(project_id: str, source_geometry_id: str, action: str, name: str = None, description: str = None, model: int = MODEL_ORB_D3_V2, extra_simulation_kwargs: dict = None):
    """
    Create a Phonon Run

    Args:
        project_id: str - The ID of the Project
        source_geometry_id: str - The ID of the Source Geometry
        action: str - The action to perform
        name: str - The name of the Phonon Run
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
    """
    Get the files associated with a Phonon Simulation Run
    """
    result = get(f"api/phonon-run-simulation-file/?phono3py_run_simulation__id={phonon_sim_run_id}")
    return result
