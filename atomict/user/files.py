from atomict.api import post


def upload_single_file(full_path: str, file_name: str, project_uuid: str = None):

    payload = {
        'users_name': file_name
    }

    if project_uuid:
        payload['project_uuid'] = project_uuid

    with open(full_path, "rb") as f:
        result = post("user/file_upload/", files={file_name: f}, payload=payload)
        return result
