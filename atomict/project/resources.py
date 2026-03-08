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


def get_project_molecule(molecule_id: str, **params) -> dict:
    return retrieve_resource("api/project-molecule", molecule_id, **params)


def list_project_molecules(**params) -> dict:
    return list_resources("api/project-molecule", **params)


def create_project_molecule(project_id: str, **fields) -> dict:
    payload = {"project": project_id}
    payload.update(fields)
    return create_resource("api/project-molecule", payload)


def update_project_molecule(
    molecule_id: str, fields: dict[str, object]
) -> dict:
    return update_resource("api/project-molecule", molecule_id, fields)


def delete_project_molecule(molecule_id: str) -> dict:
    return delete_resource("api/project-molecule", molecule_id)


def get_project_workbench_layout(layout_id: str, **params) -> dict:
    return retrieve_resource("api/project-workbench-layout", layout_id, **params)


def list_project_workbench_layouts(**params) -> dict:
    return list_resources("api/project-workbench-layout", **params)


def create_project_workbench_layout(
    project_id: str, name: str, layout: dict[str, object]
) -> dict:
    return create_resource(
        "api/project-workbench-layout",
        {"project": project_id, "name": name, "layout": layout},
    )


def update_project_workbench_layout(
    layout_id: str, fields: dict[str, object]
) -> dict:
    return update_resource("api/project-workbench-layout", layout_id, fields)


def delete_project_workbench_layout(layout_id: str) -> dict:
    return delete_resource("api/project-workbench-layout", layout_id)
