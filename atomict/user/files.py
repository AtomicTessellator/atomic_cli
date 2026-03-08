import os

from atomict.api import get, post
from atomict.resource_helpers import delete_resource, list_resources, retrieve_resource, update_resource


def upload_single_file(full_path: str, file_name: str = None, project_id: str = None):

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
    return retrieve_resource("api/user-upload", upload_id, **params)


def list_user_uploads(**params) -> dict:
    return list_resources("api/user-upload", **params)


def update_user_upload(upload_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/user-upload", upload_id, fields)


def delete_user_upload(upload_id: str) -> dict:
    return delete_resource("api/user-upload", upload_id)


def download_file(user_upload_id: str, destination_path: str):
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
    return get_user_upload(upload_id, include_content="true")
