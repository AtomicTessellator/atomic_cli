from atomict.api import get, post, delete


def get_md_simulation(simulation_id: str):
    """Get molecular dynamics simulation details.
    
    Args:
        simulation_id (str): The simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    return get(f"api/md/{simulation_id}/")


def get_md_simulation_file(file_id: str):
    """Get molecular dynamics simulation file details.
    
    Args:
        file_id (str): The file identifier.
    
    Returns:
        dict: The API response payload.
    """
    return get(f"api/md-file/{file_id}/")


def create_md_simulation(data: dict):
    """Create a new molecular dynamics simulation.
    
    Args:
        data (dict): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/md/", data)


def create_md_simulation_file(data: dict):
    """Create a new molecular dynamics simulation file.
    
    Args:
        data (dict): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/md-file/", data)


def associate_user_upload_with_md_simulation(user_upload_id: str, simulation_id: str):
    """Associate a user upload with a molecular dynamics simulation.
    
    Args:
        user_upload_id (str): The user upload identifier.
        simulation_id (str): The simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    return post(
        "api/md-file/",
        payload={"user_upload_id": user_upload_id, "md_id": simulation_id},
    )


def delete_md_simulation(simulation_id: str):
    """Delete molecular dynamics simulation.
    
    Args:
        simulation_id (str): The simulation identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete(f"api/md/{simulation_id}/")


def delete_md_simulation_file(file_id: str):
    """Delete molecular dynamics simulation file.
    
    Args:
        file_id (str): The file identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete(f"api/md-file/{file_id}/")
