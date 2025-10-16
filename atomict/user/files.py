import os
import logging

from atomict.api import get, post
from atomict.utils.path_security import (
    ensure_safe_read_path,
    ensure_safe_write_path,
    get_safe_filename,
)

logger = logging.getLogger(__name__)


def upload_single_file(full_path: str, file_name: str = None, project_id: str = None):
    """
    Upload a single file to the server.
    
    Args:
        full_path: Path to the file to upload (validated for security)
        file_name: Optional custom name for the file
        project_id: Optional project ID to associate with the file
        
    Returns:
        API response from the upload
        
    Raises:
        PathTraversalError: If the path contains traversal attempts
        FileNotFoundError: If the file doesn't exist
    """
    # Validate the file path for safe reading
    # Note: We don't restrict to a base_dir here since users should be able to
    # upload files from anywhere they have permission to read
    safe_path = ensure_safe_read_path(full_path)
    
    payload = {}

    if file_name:
        # Sanitize the filename
        payload['users_name'] = get_safe_filename(file_name)
    else:
        payload['users_name'] = os.path.basename(safe_path)

    if project_id:
        payload['project_id'] = project_id

    logger.info(f"Uploading file from: {safe_path}")
    with open(safe_path, "rb") as f:
        result = post("user/file_upload/", files={file_name: f}, payload=payload)
        return result


def download_file(user_upload_id: str, destination_path: str, base_download_dir: str = None):
    """
    Download a file from the server.
    
    Args:
        user_upload_id: ID of the file to download
        destination_path: Where to save the file
        base_download_dir: Optional base directory to restrict downloads to.
                          If not provided, uses current working directory.
        
    Returns:
        The downloaded content
        
    Raises:
        PathTraversalError: If the destination path attempts traversal
    """
    # Get content from API
    content = get(f"user/file_upload_get/{user_upload_id}/")

    # If no base directory specified, use current working directory
    if base_download_dir is None:
        base_download_dir = os.getcwd()
    
    # Validate and ensure safe write path
    # This will prevent writing outside the allowed directory
    safe_destination = ensure_safe_write_path(
        destination_path,
        base_download_dir,
        create_dirs=True
    )
    
    logger.info(f"Downloading file to: {safe_destination}")
    with open(safe_destination, "wb") as f:
        f.write(content)
    
    return content
