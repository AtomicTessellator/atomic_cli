from atomict.api import post
from atomict.resource_helpers import (
    delete_resource,
    list_resources,
    retrieve_resource,
    update_resource,
)


def get_simulation(simulation_id: str, full: bool = False, **params):
    """
    Get a SQS simulation

    Args:
        simulation_id: str - The ID of the simulation
        full: bool - Whether to get the full simulation details (default: False)
        **params: Additional GET parameters to pass to the API
    """
    query_params = params.copy()
    if full:
        query_params["full"] = "true"
    return retrieve_resource("api/sqs-exploration", simulation_id, **query_params)


def get_sqs_exploration(exploration_id: str, **params):
    return get_simulation(exploration_id, **params)


def list_sqs_explorations(**params):
    return list_resources("api/sqs-exploration", **params)


def create_sqs_exploration(payload: dict[str, object]):
    return post("api/sqs-exploration/", payload=payload)


def update_sqs_exploration(exploration_id: str, fields: dict[str, object]):
    return update_resource("api/sqs-exploration", exploration_id, fields)


def delete_sqs_exploration(exploration_id: str):
    return delete_resource("api/sqs-exploration", exploration_id)


def associate_user_upload_with_sqs_simulation(user_upload_id: str, exploration_id: str):
    """
    Associate a user upload with a SQS simulation
    """
    result = post(
        "api/sqs-simulation-file/",
        payload={"user_upload_id": user_upload_id, "exploration_id": exploration_id},
    )
    return result


def get_sqs_target_concentration(target_concentration_id: str, **params):
    return retrieve_resource(
        "api/sqs-target-concentration", target_concentration_id, **params
    )


def list_sqs_target_concentrations(**params):
    return list_resources("api/sqs-target-concentration", **params)


def create_sqs_target_concentration(payload: dict[str, object]):
    return post("api/sqs-target-concentration/", payload=payload)


def update_sqs_target_concentration(
    target_concentration_id: str, fields: dict[str, object]
):
    return update_resource(
        "api/sqs-target-concentration", target_concentration_id, fields
    )


def delete_sqs_target_concentration(target_concentration_id: str):
    return delete_resource("api/sqs-target-concentration", target_concentration_id)


def get_sqs_simulation_file(file_id: str, **params):
    return retrieve_resource("api/sqs-simulation-file", file_id, **params)


def list_sqs_simulation_files(**params):
    return list_resources("api/sqs-simulation-file", **params)


def create_sqs_simulation_file(payload: dict[str, object]):
    return post("api/sqs-simulation-file/", payload=payload)


def update_sqs_simulation_file(file_id: str, fields: dict[str, object]):
    return update_resource("api/sqs-simulation-file", file_id, fields)


def delete_sqs_simulation_file(file_id: str):
    return delete_resource("api/sqs-simulation-file", file_id)
