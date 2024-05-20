from atomict.api import post


def upload_single_file(full_path: str, file_name: str):
    with open(full_path, "rb") as f:
        result = post(
            "user/file_upload/",
            files={file_name: f},
            payload={'users_name': file_name}
        )
        return result
