from atomict.resource_helpers import (
    create_resource,
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_organisation(organisation_id: str, **params) -> dict:
    """Get organisation details.
    
    Args:
        organisation_id (str): The organisation identifier.
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return retrieve_resource("api/organisation", organisation_id, **params)


def list_organisations(**params) -> dict:
    """List organisations.
    
    Args:
        **params (Any): Additional query parameters to include in the request.
    
    Returns:
        dict: The API response payload.
    """
    return list_resources("api/organisation", **params)


def create_organisation(payload: dict[str, object]) -> dict:
    """Create a new organisation.
    
    Args:
        payload (dict[str, object]): The payload to send to the API.
    
    Returns:
        dict: The API response payload.
    """
    return create_resource("api/organisation", payload)


def update_organisation(organisation_id: str, fields: dict[str, object]) -> dict:
    """Update organisation.
    
    Args:
        organisation_id (str): The organisation identifier.
        fields (dict[str, object]): Field values to update on the resource.
    
    Returns:
        dict: The API response payload.
    """
    return update_resource("api/organisation", organisation_id, fields)


def delete_organisation(organisation_id: str) -> dict:
    """Delete organisation.
    
    Args:
        organisation_id (str): The organisation identifier.
    
    Returns:
        dict: The API response payload.
    """
    return delete_resource("api/organisation", organisation_id)


def get_organisation_invite(invite_id: str, **params) -> dict:
    return retrieve_resource("api/organisation-invite", invite_id, **params)


def list_organisation_invites(**params) -> dict:
    return list_resources("api/organisation-invite", **params)


def create_organisation_invite(payload: dict[str, object]) -> dict:
    return create_resource("api/organisation-invite", payload)


def update_organisation_invite(invite_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/organisation-invite", invite_id, fields)


def delete_organisation_invite(invite_id: str) -> dict:
    return delete_resource("api/organisation-invite", invite_id)


def get_organisation_user(user_id: str, **params) -> dict:
    return retrieve_resource("api/organisation-user", user_id, **params)


def list_organisation_users(**params) -> dict:
    return list_resources("api/organisation-user", **params)


def create_organisation_user(payload: dict[str, object]) -> dict:
    return create_resource("api/organisation-user", payload)


def update_organisation_user(user_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/organisation-user", user_id, fields)


def delete_organisation_user(user_id: str) -> dict:
    return delete_resource("api/organisation-user", user_id)
