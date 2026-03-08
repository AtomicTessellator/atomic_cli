from atomict.api import get
from atomict.resource_helpers import delete_resource, list_resources, retrieve_resource, update_resource


def get_user(user_id: str, **params) -> dict:
    return retrieve_resource("api/user", user_id, **params)


def list_users(**params) -> dict:
    return list_resources("api/user", **params)


def update_user(user_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/user", user_id, fields)


def delete_user(user_id: str) -> dict:
    return delete_resource("api/user", user_id)


def lookup_simulation_workspace(simulation_uuid: str, **params) -> dict:
    query = dict(params)
    return get(f"simulation/lookup/{simulation_uuid}/", params=query or None)
