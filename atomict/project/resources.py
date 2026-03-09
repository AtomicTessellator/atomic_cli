from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_project_star(star_id: str, **params) -> dict:
    return retrieve_resource("api/project-star", star_id, **params)


def list_project_stars(**params) -> dict:
    return list_resources("api/project-star", **params)


def create_project_star(project_id: str) -> dict:
    return create_resource("api/project-star", {"project": project_id})


def update_project_star(star_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/project-star", star_id, fields)


def delete_project_star(star_id: str) -> dict:
    return delete_resource("api/project-star", star_id)


def get_project_workbench_layout(layout_id: str, **params) -> dict:
    """Get project workbench layout details.
    
    Args:
        layout_id (str): The layout identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/project-workbench-layout", layout_id, **params)


def list_project_workbench_layouts(**params) -> dict:
    """List project workbench layouts.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/project-workbench-layout", **params)


def create_project_workbench_layout(
    project_id: str, name: str, layout: dict[str, object]
) -> dict:
    """Create a new project workbench layout.
    
    Args:
        project_id (str): The project identifier.
        name (str): The resource name.
        layout (dict[str, object]): The layout.
    
    Returns:
        dict: The API response payload.
    """
    return create_resource(
        "api/project-workbench-layout",
        {"project": project_id, "name": name, "layout": layout},
    )


def update_project_workbench_layout(
    layout_id: str, fields: dict[str, object]
) -> dict:
    """Update project workbench layout.
    
    Args:
        layout_id (str): The layout identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/project-workbench-layout", layout_id, fields)


def delete_project_workbench_layout(layout_id: str) -> dict:
    """Delete project workbench layout.
    
    Args:
        layout_id (str): The layout identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/project-workbench-layout", layout_id)
