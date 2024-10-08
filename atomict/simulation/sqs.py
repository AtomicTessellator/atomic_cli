from atomict.api import get, post


def get_simulation(simulation_id: str, full: bool = False):
    """
    Get a SQS simulation
    """
    if full:
        result = get(f"api/sqs-exploration/{simulation_id}/?full=true")
    else:
        result = get(f"api/sqs-exploration/{simulation_id}/")

    return result


def associate_user_upload_with_sqs_simulation(
    user_upload_id: str, sqs_simulation_id: str
):
    """
    Associate a user upload with a SQS simulation
    """
    result = post(
        "api/sqs-simulation-file/",
        payload={"user_upload": user_upload_id, "simulation": sqs_simulation_id},
    )