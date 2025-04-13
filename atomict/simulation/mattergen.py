from atomict.api import get, post


def get_mattergen(id: str):
    """
    Get MatterGen
    """
    result = get(f"api/mattergen-exploration/{id}/")
    return result


def associate_user_upload_with_mattergen(user_upload_id: str, mattergen_id: str):
    """
    Associate a user upload with a MatterGen
    """
    result = post(
        "api/mattergen-exploration-file/",
        payload={"user_upload_id": user_upload_id, "mattergen_id": mattergen_id},
    )
    return result
