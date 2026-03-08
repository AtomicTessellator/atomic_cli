from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def create_project_note(
    project_id: str, title: str, content: str, show_description: bool = True
) -> dict:

    payload = {
        "project": project_id,
        "title": title,
        "content_html": content,
        "show_description": show_description,
    }

    return create_resource("api/project-note", payload)


def get_project_note(note_id: str, **params) -> dict:
    return retrieve_resource("api/project-note", note_id, **params)


def list_project_notes(**params) -> dict:
    return list_resources("api/project-note", **params)


def update_project_note(note_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/project-note", note_id, fields)


def delete_project_note(note_id: str) -> dict:
    return delete_resource("api/project-note", note_id)
