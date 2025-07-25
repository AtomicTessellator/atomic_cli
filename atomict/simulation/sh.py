from atomict.api import get, post


def get_sh(id: str, **params):
    """
    Get SH
    """
    result = get(f"api/encoding-task/{id}/", params=params)
    return result


def associate_user_upload_with_sh(user_upload_id: str, sh_id: str):
    """
    Associate a user upload with a SH
    """
    result = post(
        "api/encoding-file/",
        payload={"user_upload_id": user_upload_id, "encoding_id": sh_id},
    )
    return result


def get_sh_files(sh_id: str):
    """
    Get the files associated with a SH
    """
    result = get(f"api/encoding-file/?encoding__id={sh_id}")
    return result
