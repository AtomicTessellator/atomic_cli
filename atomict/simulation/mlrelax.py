from atomict.api import get, post, patch


def get_mlrelax(id: str):
    """
    Get MLRelaxation
    """
    result = get(f"api/mlrelax/{id}/")
    return result


def associate_user_upload_with_mlrelaxation(user_upload_id: str, mlrelax_id: str):
    """
    Associate a user upload with a MLRelaxation
    """
    result = post(
        "api/mlrelax-file/",
        payload={"user_upload": user_upload_id, "mlrelax": mlrelax_id},
    )
    return result
