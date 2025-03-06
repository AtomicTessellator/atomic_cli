from atomict.api import get, post


def get_defect_exploration(exploration_id: str):
    return get(f"api/defect-exploration/{exploration_id}/")


def get_defect_exploration_file(file_id: str):
    return get(f"api/defect-exploration-file/{file_id}/")


def create_defect_exploration(data: dict):
    return post("api/defect-exploration/", data)


def create_defect_exploration_file(data: dict):
    return post("api/defect-exploration-file/", data)


def associate_user_upload_with_defect_exploration(user_upload_id: str, exploration_id: str):
    return post(
        "api/defect-exploration-file/",
        payload={"user_upload_id": user_upload_id, "analysis_id": exploration_id},
    )
