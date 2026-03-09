import os

from atomict.api import get, post
from atomict.resource_helpers import delete_resource, list_resources, retrieve_resource, update_resource


def upload_single_file(full_path: str, file_name: str = None, project_id: str = None):

    """Upload a user file.
    
    Args:
        full_path (str): The path to the local file to upload.
        file_name (str | None): The display name to use for the uploaded file.
        project_id (str | None): The project identifier.
    
    Returns:
        dict: The API response payload.
    """
    payload = {}

    if file_name:
        payload['users_name'] = file_name
    else:
        payload['users_name'] = os.path.basename(full_path)

    if project_id:
        payload['project_id'] = project_id

    with open(full_path, "rb") as f:
        result = post("user/file_upload/", files={file_name: f}, payload=payload)
        return result


def get_user_upload(upload_id: str, **params) -> dict:
    """Get user upload details.
    
    Args:
        upload_id (str): The user upload identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/user-upload", upload_id, **params)


def list_user_uploads(**params) -> dict:
    """List user uploads.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/user-upload", **params)


def update_user_upload(upload_id: str, fields: dict[str, object]) -> dict:
    """Update user upload.
    
    Args:
        upload_id (str): The user upload identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/user-upload", upload_id, fields)


def delete_user_upload(upload_id: str) -> dict:
    """Delete user upload.
    
    Args:
        upload_id (str): The user upload identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/user-upload", upload_id)


def download_file(user_upload_id: str, destination_path: str):
    """Download a user upload to a local path.
    
    Args:
        user_upload_id (str): The user upload identifier.
        destination_path (str): The local path where the downloaded file should be written.
    
    Returns:
        bytes: The downloaded binary file content.
    """
    content = get(f"user/file_upload_get/{user_upload_id}/")

    # Write the content to the destination path
    # if there's a directory path in the destination path, create the directory
    destination_dir = os.path.dirname(destination_path)
    if destination_dir:
        os.makedirs(destination_dir, exist_ok=True)

    with open(destination_path, "wb") as f:
        f.write(content)
    return content


def download_user_upload_content(upload_id: str) -> dict:
    """Download the content of a user upload.
    
    Args:
        upload_id (str): The user upload identifier.
    
    Returns:
        dict: The API response payload.
    """
    return get_user_upload(upload_id, include_content="true")
