import os
from typing import Any, Dict, List, Optional, Union

from atomict.api import delete, get, patch, post


def upload_single_file(
    full_path: str, file_name: str, project_uuid: Optional[str] = None
):

    payload = {"users_name": file_name}

    if project_uuid:
        payload["project_uuid"] = project_uuid

    with open(full_path, "rb") as f:
        result = post("user/file_upload/", files={file_name: f}, payload=payload)
        return result


def download_file(user_upload_uuid: str, destination_path: str):
    content = get(f"user/file_upload_get/{user_upload_uuid}/")

    # Write the content to the destination path
    # if there's a directory path in the destination path, create the directory
    destination_dir = os.path.dirname(destination_path)
    if destination_dir:
        os.makedirs(destination_dir, exist_ok=True)

    if content is not None:
        with open(destination_path, "wb") as f:
            f.write(content)
    return content


def delete_user_upload(upload_uuid: str) -> bool:
    """Delete a user upload by UUID.

    Args:
        upload_uuid: UUID of the upload to delete

    Returns:
        True if deletion was successful

    Raises:
        APIValidationError: If the upload doesn't exist or user lacks permission
    """
    result = delete(f"user/file_upload_delete/{upload_uuid}/")
    return result.get("success", True) if result is not None else False


def get_user_uploads(
    user_id: Optional[str] = None,
    search: Optional[str] = None,
    ordering: Optional[str] = None,
    limit: Optional[int] = None,
    offset: Optional[int] = None,
) -> Any:
    """List user uploads with filtering options.

    Args:
        user_id: Filter by specific user ID
        search: Search term for file names/descriptions
        ordering: Field to order results by (e.g., "-uploaded", "orig_name")
        limit: Maximum number of results to return
        offset: Number of results to skip

    Returns:
        Dictionary containing upload list and pagination info
    """
    # Build query string manually since API doesn't support params dict
    query_parts = []
    if user_id:
        query_parts.append(f"user={user_id}")
    if search:
        query_parts.append(f"search={search}")
    if ordering:
        query_parts.append(f"ordering={ordering}")
    if limit:
        query_parts.append(f"limit={limit}")
    if offset:
        query_parts.append(f"offset={offset}")

    path = "api/user-upload/"
    if query_parts:
        path += "?" + "&".join(query_parts)

    result = get(path)
    return result if result is not None else {}


def get_user_upload(upload_uuid: str, include_content: bool = False) -> Any:
    """Get user upload details, optionally with file content.

    Args:
        upload_uuid: UUID of the upload to retrieve
        include_content: Whether to include base64-encoded file content

    Returns:
        Dictionary containing upload details and optionally file content
    """
    path = f"api/user-upload/{upload_uuid}/"
    if include_content:
        path += "?include_content=true"

    result = get(path)
    return result if result is not None else {}


def update_user_upload(upload_uuid: str, **updates) -> Any:
    """Update user upload metadata.

    Args:
        upload_uuid: UUID of the upload to update
        **updates: Fields to update (e.g., users_description, users_name)

    Returns:
        Updated upload details
    """
    return patch(f"api/user-upload/{upload_uuid}/", payload=updates)


def upload_multiple_files(
    file_paths: List[str], project_uuid: Optional[str] = None
) -> List[Dict]:
    """Bulk upload multiple files.

    Args:
        file_paths: List of file paths to upload
        project_uuid: Optional project to associate uploads with

    Returns:
        List of upload results
    """
    results = []
    for file_path in file_paths:
        file_name = os.path.basename(file_path)
        try:
            result = upload_single_file(file_path, file_name, project_uuid)
            results.append({"file": file_path, "success": True, "result": result})
        except Exception as e:
            results.append({"file": file_path, "success": False, "error": str(e)})

    return results


def get_file_upload_filter(
    file_types: Optional[List[str]] = None,
    project_uuid: Optional[str] = None,
    **filters,
) -> Any:
    """Get filtered list of file uploads with advanced filtering.

    Args:
        file_types: List of file types to filter by (e.g., ["xyz", "cif"])
        project_uuid: Filter by project UUID
        **filters: Additional filter parameters

    Returns:
        Filtered upload results
    """
    params = {}
    if file_types:
        params["file_types"] = ",".join(file_types)
    if project_uuid:
        params["project_uuid"] = project_uuid
    params.update(filters)

    return post("user/file_upload_filter/", payload=params)
