from atomict.api import get, post


def upload_single_file(full_path: str, file_name: str, project_uuid: str = None):

    payload = {
        'users_name': file_name
    }

    if project_uuid:
        payload['project_uuid'] = project_uuid

    with open(full_path, "rb") as f:
        result = post("user/file_upload/", files={file_name: f}, payload=payload)
        return result


def download_file(user_upload_uuid: str, destination_path: str):
    content = get(f"user/file_upload_get/{user_upload_uuid}/")
    with open(destination_path, "wb") as f:
        f.write(content)
    return content
