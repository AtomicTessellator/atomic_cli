from atomict.api import get
from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)

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

    return create_resource("api/project-tag", {"tag": name, "color": tag_color})


def get_tag_by_name(tag: str) -> dict:
    response = get(f"api/project-tag/?tag={tag}")

    return response['results'][0]


def tag_exists(tag: str) -> bool:
    response = get(f"api/project-tag/?tag={tag}")

    return response['count'] > 0


def create_project_tag(project_id: str, tag_id: str) -> dict:
    return create_resource(
        "api/project-tag-project",
        {"project": project_id, "project_tag": tag_id},
    )


def project_tag_exists(project_id: str, tag_id: str) -> bool:
    response = get(f"api/project-tag-project/?project={project_id}&project_tag={tag_id}")

    return response['count'] > 0


def get_project_tag(tag_id: str, **params) -> dict:
    return retrieve_resource("api/project-tag", tag_id, **params)


def list_project_tags(**params) -> dict:
    return list_resources("api/project-tag", **params)


def update_project_tag(tag_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/project-tag", tag_id, fields)


def delete_project_tag(tag_id: str) -> dict:
    return delete_resource("api/project-tag", tag_id)


def get_project_tag_project(assignment_id: str, **params) -> dict:
    return retrieve_resource("api/project-tag-project", assignment_id, **params)


def list_project_tag_projects(**params) -> dict:
    return list_resources("api/project-tag-project", **params)


def update_project_tag_project(
    assignment_id: str, fields: dict[str, object]
) -> dict:
    return update_resource("api/project-tag-project", assignment_id, fields)


def delete_project_tag_project(assignment_id: str) -> dict:
    return delete_resource("api/project-tag-project", assignment_id)
