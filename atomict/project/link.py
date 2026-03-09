from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def create_object_link(
    project_id: str, src_object_id: str, dst_object_id: str) -> dict:
    """Create a new object link.

    An object link is a relationship between two objects in a project.
    It is used to connect objects together, such as a note to a file or a simulation to a note.
    This is used to layout visually the project structure, it does not affect the data model.
    
    Args:
        project_id (str): The project UUID.
        src_object_id (str): The source object UUID.
        dst_object_id (str): The destination object UUID.
    
    Returns:
        dict: The API response payload.
    """
    payload = {
        "project": project_id,
        "src_id": src_object_id,
        "dst_id": dst_object_id
    }

    return create_resource("api/object-link", payload)


def delete_object_link(link_id: str) -> dict:
    """Delete object link.
    
    Args:
        link_id (str): The object link UUID.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/object-link", link_id)


def get_object_link(link_id: str, **params) -> dict:
    """Get object link details.
    
    Args:
        link_id (str): The object link UUID.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/object-link", link_id, **params)


def list_object_links(**params) -> dict:
    """List object links.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/object-link", **params)


def update_object_link(link_id: str, fields: dict[str, object]) -> dict:
    """Update object link.
    
    Args:
        link_id (str): The object link UUID.
        fields (dict[str, object]): Field values to update on the resource.

        The only fields that can be updated are:
        - src_id
        - dst_id
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/object-link", link_id, fields)
