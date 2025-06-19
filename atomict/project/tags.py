from atomict.api import get, post, patch, delete

VALID_TAG_COLOURS = [
    "bg-success",
    "bg-primary",
    "bg-secondary",
    "bg-danger",
    "bg-warning",
    "bg-info",
    "bg-light",
    "bg-dark"
]

def create_tag(name: str, tag_color: str) -> int:

    if tag_color not in VALID_TAG_COLOURS:
        raise ValueError(f"Invalid tag color: {tag_color}, choose from {VALID_TAG_COLOURS}")

    response = post(
        "api/project-tag/", {"tag": name, "color": tag_color}, extra_headers={"Content-Type": "application/json"})
    return response


def get_tag_by_name(tag: str) -> dict:
    response = get(f"api/project-tag/?tag={tag}")

    return response['results'][0]


def tag_exists(tag: str) -> bool:
    response = get(f"api/project-tag/?tag={tag}")

    return response['count'] > 0


def create_project_tag(project_id: str, tag_id: str) -> dict:
    response = post(
        "api/project-tag-project/",
        {"project": project_id, "project_tag": tag_id},
        extra_headers={"Content-Type": "application/json"},
    )
    return response


def project_tag_exists(project_id: str, tag_id: str) -> bool:
    response = get(f"api/project-tag-project/?project={project_id}&project_tag={tag_id}")

    return response['count'] > 0


def delete_project_tag(tag_id: str) -> dict:
    response = delete(f"api/project-tag/{tag_id}/")
    return response


def list_project_tags() -> dict:
    response = get("api/project-tag/")
    return response


def get_project_tag(tag_id: str) -> dict:
    response = get(f"api/project-tag/{tag_id}/")
    return response


def update_project_tag(tag_id: str, tag: str = None, color: str = None) -> dict:
    payload = {}
    
    if tag is not None:
        payload["tag"] = tag
    
    if color is not None:
        if color not in VALID_TAG_COLOURS:
            raise ValueError(f"Invalid tag color: {color}, choose from {VALID_TAG_COLOURS}")
        payload["color"] = color
    
    if not payload:
        raise ValueError("At least one parameter (tag or color) must be provided")
    
    response = patch(f"api/project-tag/{tag_id}/", payload)
    return response


def delete_project_tag_association(association_id: str) -> dict:
    response = delete(f"api/project-tag-project/{association_id}/")
    return response


def list_project_tag_associations(project_id: str = None, tag_id: str = None) -> dict:
    params = []
    
    if project_id is not None:
        params.append(f"project={project_id}")
    
    if tag_id is not None:
        params.append(f"project_tag={tag_id}")
    
    query_string = "?" + "&".join(params) if params else ""
    response = get(f"api/project-tag-project/{query_string}")
    return response
