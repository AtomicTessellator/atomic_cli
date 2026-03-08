from atomict.api import delete, get, post
from atomict.resource_helpers import list_resources, retrieve_resource, update_resource


def create_project(name: str, description: str = None) -> dict:

    payload = {
        "name": name,
        "description_html": description,
    }

    response = post(
        "api/project/", payload, extra_headers={"Content-Type": "application/json"})
    return response


def get_project(project_id: str, **params) -> dict:
    return retrieve_resource("api/project", project_id, **params)


def list_projects(**params) -> dict:
    return list_resources("api/project", **params)


def update_project(project_id: str, fields: dict[str, object]) -> dict:
    return update_resource("api/project", project_id, fields)


def delete_project(project_id: str) -> dict:
    response = delete(f"api/project/{project_id}/")

    return response


def project_exists(name: str) -> bool:
    response = get(f"api/project/?name={name}")

    return response['count'] > 0


def get_project_by_name(name: str) -> dict:
    response = get(f"api/project/?name={name}")

    return response['results'][0]
