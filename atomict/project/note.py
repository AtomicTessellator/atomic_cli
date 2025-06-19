from typing import Optional
from atomict.api import delete, get, patch, post


def create_project_note(
    project_id: str, title: str, content: str, show_description: bool = True
):
    """
    Create a project note.
    
    Note: Backend automatically sets show_description=True when content is not empty,
    so show_description=False only works for notes with empty content.
    
    Args:
        project_id: The UUID of the project
        title: The title of the note
        content: The HTML content of the note
        show_description: Whether to show description (limited by backend logic)
        
    Returns:
        dict: Created note data from API
    """
    payload = {
        "project": project_id,
        "title": title,
        "content_html": content,
        "show_description": show_description,
    }

    response = post(
        "api/project-note/", payload, extra_headers={"Content-Type": "application/json"}
    )
    return response


def delete_project_note(note_id: str):
    response = delete(f"api/project-note/{note_id}/")
    return response


def get_project_note(note_id: str):
    response = get(f"api/project-note/{note_id}/")
    return response


def list_project_notes(project_id: Optional[str] = None):
    url = "api/project-note/"
    if project_id:
        url += f"?project__id={project_id}"  # Backend uses project__id filter
    response = get(url)
    return response


def update_project_note(
    note_id: str, 
    title: Optional[str] = None, 
    content: Optional[str] = None, 
    show_description: Optional[bool] = None
):
    """
    Update a project note.
    
    Note: Backend automatically sets show_description=True when content is not empty,
    so show_description=False only works for notes with empty content.
    
    Args:
        note_id: The UUID of the note to update
        title: Optional new title for the note
        content: Optional new HTML content for the note
        show_description: Optional flag for whether to show description (limited by backend logic)
        
    Returns:
        dict: Updated note data from API
        
    Raises:
        APIValidationError: If note_id is invalid or validation fails
    """
    payload = {}
    
    # Backend requires either title or content to be non-empty
    # If updating only show_description, we need to include current values
    if title is None and content is None and show_description is not None:
        # Get current note to preserve title/content during show_description update
        current_note = get_project_note(note_id)
        if isinstance(current_note, dict):
            if current_note.get("title"):
                payload["title"] = current_note["title"]
            if current_note.get("content_html"):
                payload["content_html"] = current_note["content_html"]
    
    if title is not None:
        payload["title"] = title
    if content is not None:
        payload["content_html"] = content
    if show_description is not None:
        payload["show_description"] = show_description

    response = patch(f"api/project-note/{note_id}/", payload)
    return response
