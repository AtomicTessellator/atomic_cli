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

    """Create a new project note.
    
    Args:
        project_id (str): The project identifier.
        title (str): The note title.
        content (str): The note content, HTML is supported.
        show_description (bool): Whether the description should be shown with the note.
    
    Returns:
        dict: The API response payload.
    """
    payload = {
        "project": project_id,
        "title": title,
        "content_html": content,
        "show_description": show_description,
    }

    return create_resource("api/project-note", payload)


def get_project_note(note_id: str, **params) -> dict:
    """Get project note details.
    
    Args:
        note_id (str): The note identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/project-note", note_id, **params)


def list_project_notes(**params) -> dict:
    """List project notes.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/project-note", **params)


def update_project_note(note_id: str, fields: dict[str, object]) -> dict:
    """Update project note.
    
    Args:
        note_id (str): The note identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/project-note", note_id, fields)


def delete_project_note(note_id: str) -> dict:
    """Delete project note.
    
    Args:
        note_id (str): The note identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/project-note", note_id)
