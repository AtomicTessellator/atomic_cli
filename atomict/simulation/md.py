from atomict.api import get, post, delete


def get_md_simulation(simulation_id: str):
    return get(f"api/md/{simulation_id}/")


def get_md_simulation_file(file_id: str):
    return get(f"api/md-file/{file_id}/")


def create_md_simulation(data: dict):
    return post("api/md/", data)


def create_md_simulation_file(data: dict):
    return post("api/md-file/", data)


def associate_user_upload_with_md_simulation(user_upload_id: str, simulation_id: str):
    return post(
        "api/md-file/",
        payload={"user_upload_id": user_upload_id, "md_id": simulation_id},
    )


def delete_md_simulation(simulation_id: str):
    return delete(f"api/md/{simulation_id}/")


def delete_md_simulation_file(file_id: str):
    return delete(f"api/md-file/{file_id}/")
