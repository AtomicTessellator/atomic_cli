from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def create_object_link(
    project_id: str, src_object_id: str, dst_object_id: str) -> dict:
    """
    Create a link between an object and a project
    """
    payload = {
        "project": project_id,
        "src_id": src_object_id,
        "dst_id": dst_object_id
    }

    return create_resource("api/object-link", payload)


def delete_object_link(link_id: str) -> dict:
    """
    Delete a link between an object and a project
    """
    return delete_resource("api/object-link", link_id)


def get_object_link(link_id: str, **params) -> dict:
    return retrieve_resource("api/object-link", link_id, **params)


def list_object_links(**params) -> dict:
    return list_resources("api/object-link", **params)


def update_object_link(link_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/object-link", link_id, fields)
