from atomict.api import post, delete


def create_object_link(
    project_id: str, src_object_id: str, dst_object_id: str) -> dict:
    """
    Create a link between an object and a project
    """
    payload = {
        "project": project_id,
        "src_id": src_object_id,
        "dst_id": dst_object_id
    }

    response = post(
        "api/object-link/", payload, extra_headers={"Content-Type": "application/json"}
    )
    return response


def delete_object_link(link_id: str) -> dict:
    """
    Delete a link between an object and a project
    """
    response = delete(f"api/object-link/{link_id}/")
    return response
