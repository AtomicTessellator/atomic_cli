from atomict.api import get, post


def get_sh(id: str, **params):
    """Get encoding task details.
    
    Args:
        id (str): The resource identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    result = get(f"api/encoding-task/{id}/", params=params)
    return result


def associate_user_upload_with_sh(user_upload_id: str, sh_id: str):
    """Associate a user upload with an encoding task.
    
    Args:
        user_upload_id (str): The user upload identifier.
        sh_id (str): The encoding task identifier.
    
    Returns:
        dict: The API response payload.
    """
    result = post(
        "api/encoding-file/",
        payload={"user_upload_id": user_upload_id, "encoding_id": sh_id},
    )
    return result


def get_sh_files(sh_id: str):
    """List files associated with an encoding task.
    
    Args:
        sh_id (str): The encoding task identifier.
    
    Returns:
        dict: The API response payload.
    """
    result = get(f"api/encoding-file/?encoding__id={sh_id}")
    return result
