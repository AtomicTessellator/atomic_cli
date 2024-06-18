from atomict.api import get, post

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
