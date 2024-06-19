from atomict.api import post


def create_project_note(
    project_id: str, title: str, content: str, show_description: bool = True
) -> dict:

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
