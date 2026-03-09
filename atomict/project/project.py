from atomict.api import delete, get, post
from atomict.resource_helpers import list_resources, retrieve_resource, update_resource


def create_project(name: str, description: str = None) -> dict:

    """Create a new project.

    A project is a collection of resources, such as notes, files (e.g. CIF files)
    and simulations.
    
    Args:
        name (str): The short, descriptive name of the project, e.g. "Mechanical Properties of Silicon".
        description (str | None): The project description, 5 lines max.
    
    Returns:
        dict: The API response payload.
    """
    payload = {
        "name": name,
        "description_html": description,
    }

    response = post(
        "api/project/", payload, extra_headers={"Content-Type": "application/json"})
    return response


def get_project(project_id: str, **params) -> dict:
    """Get project details.
    
    Args:
        project_id (str): The project UUID.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/project", project_id, **params)


def list_projects(**params) -> dict:
    """List projects.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/project", **params)


def update_project(project_id: str, fields: dict[str, object]) -> dict:
    """Update project.
    
    Args:
        project_id (str): The project UUID.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/project", project_id, fields)


def delete_project(project_id: str) -> dict:
    """Delete project.
    
    Args:
        project_id (str): The project UUID.
    
    Returns:
        dict: The API response payload.
    """
    response = delete(f"api/project/{project_id}/")

    return response


def project_exists(name: str) -> bool:
    """Check whether a project exists.
    
    Args:
        name (str): The project name.
    
    Returns:
        bool: True if the requested condition is met; otherwise, False.
    """
    response = get(f"api/project/?name={name}")

    return response['count'] > 0


def get_project_by_name(name: str) -> dict:
    """Get a project by name.
    
    Args:
        name (str): The project name.
    
    Returns:
        dict: The API response payload.
    """
    response = get(f"api/project/?name={name}")

    return response['results'][0]
