from atomict.api import get, post


def get_defect_exploration(exploration_id: str):
    """Get defect exploration details.
    
    Args:
        exploration_id (str): The exploration identifier.
    
    Returns:
        dict: The API response payload.
    """
    return get(f"api/defect-exploration/{exploration_id}/")


def get_defect_exploration_file(file_id: str):
    """Get defect exploration file details.
    
    Args:
        file_id (str): The file identifier.
    
    Returns:
        dict: The API response payload.
    """
    return get(f"api/defect-exploration-file/{file_id}/")


def create_defect_exploration(data: dict):
    """Create a new defect exploration.
    
    Args:
        data (dict): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/defect-exploration/", data)


def create_defect_exploration_file(data: dict):
    """Create a new defect exploration file.
    
    Args:
        data (dict): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return post("api/defect-exploration-file/", data)


def associate_user_upload_with_defect_exploration(user_upload_id: str, exploration_id: str):
    """Associate a user upload with a defect exploration.
    
    Args:
        user_upload_id (str): The user upload identifier.
        exploration_id (str): The exploration identifier.
    
    Returns:
        dict: The API response payload.
    """
    return post(
        "api/defect-exploration-file/",
        payload={"user_upload_id": user_upload_id, "analysis_id": exploration_id},
    )
