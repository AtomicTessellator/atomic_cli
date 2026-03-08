from typing import Union

from atomict.api import post
from atomict.infra.distwork.task import SimulationAction
from atomict.resource_helpers import (
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)
from atomict.simulation.models import MODEL_ORB_V3_CONSERVATIVE, MODEL_ESEN_30M_OAM, MODEL_UMA

COMPUTATION_TYPE_RELAXATION = 0
COMPUTATION_TYPE_SINGLE_POINT = 1


def get_mlrelax(id: str, **params):
    """
    Get MLRelaxation

    Args:
        id: str - The ID of the MLRelaxation
        **params: Additional GET parameters to pass to the API
    """
    return retrieve_resource("api/mlrelax", id, **params)


def get_mlrelaxation(mlrelax_id: str, **params):
    return get_mlrelax(mlrelax_id, **params)


def list_mlrelaxations(**params):
    return list_resources("api/mlrelax", **params)


def associate_user_upload_with_mlrelaxation(user_upload_id: str, mlrelax_id: str):
    """
    Associate a user upload with a MLRelaxation
    """
    result = post(
        "api/mlrelax-file/",
        payload={"user_upload_id": user_upload_id, "mlrelax_id": mlrelax_id},
    )
    return result


def create_mlrelaxation(
    project_id: str,
    source_geometry_id: str,
    action: SimulationAction,
    name: str = None,
    description: str = None,
    computation_type: int = COMPUTATION_TYPE_RELAXATION,
    f_max: float = None,
    model: int = MODEL_ESEN_30M_OAM,
    calculator: Union[int, None] = None,
    extra_simulation_kwargs: dict = None,
):
    """
    Create a MLRelaxation
    """

    if action not in [SimulationAction.SAVE_DRAFT, SimulationAction.LAUNCH]:
        raise ValueError("Action must be 'SimulationAction.SAVE_DRAFT' or 'SimulationAction.LAUNCH'")

    if computation_type not in [
        COMPUTATION_TYPE_RELAXATION,
        COMPUTATION_TYPE_SINGLE_POINT,
    ]:
        raise ValueError(
            "Invalid computation type. Please use atomict.simulation.mlrelax.COMPUTATION_TYPE_RELAXATION or atomict.simulation.mlrelax.COMPUTATION_TYPE_SINGLE_POINT."
        )

    payload = {
        "project": project_id,
        "source_geometry_id": source_geometry_id,
        "action": action,
        "name": name,
        "description": description,
        "computation_type": computation_type,
        "f_max": f_max,
    }

    # backwards compat - accept 'model' param but always send 'calculator' to API
    # (model field was renamed to calculator in the backend)
    if calculator is not None:
        payload["calculator"] = calculator
    else:
        payload["calculator"] = model

    if extra_simulation_kwargs:
        payload.update(extra_simulation_kwargs)

    result = post("api/mlrelax/", payload)
    return result


def update_mlrelaxation(mlrelax_id: str, fields: dict[str, object]):
    return update_resource("api/mlrelax", mlrelax_id, fields)


def delete_mlrelaxation(mlrelax_id: str):
    return delete_resource("api/mlrelax", mlrelax_id)


def get_mlrelax_files(mlrelax_id: str):
    """
    Get the files associated with a MLRelaxation
    """
    return list_resources("api/mlrelax-file", mlrelax__id=mlrelax_id)


def get_mlrelax_file(file_id: str, **params):
    return retrieve_resource("api/mlrelax-file", file_id, **params)


def list_mlrelax_files(**params):
    return list_resources("api/mlrelax-file", **params)


def create_mlrelax_file(payload: dict[str, object]):
    return post("api/mlrelax-file/", payload=payload)


def update_mlrelax_file(file_id: str, fields: dict[str, object]):
    return update_resource("api/mlrelax-file", file_id, fields)


def delete_mlrelax_file(file_id: str):
    return delete_resource("api/mlrelax-file", file_id)
