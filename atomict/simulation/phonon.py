from atomict.api import post


def associate_user_upload_with_phonon_analysis(
    user_upload_id: str, pa_id: str
):
    """
    Associate a user upload with a Phonon Analysis
    """
    result = post(
        "api/phonon-analysis-file/",
        payload={"user_upload": user_upload_id, "analysis": pa_id},
    )
    return result
